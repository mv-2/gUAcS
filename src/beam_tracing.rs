use crate::interface::{BeamConfigRust, EnvConfig, IsoSpace, ProgConfig, SolverMethod};
use crate::math_util::Mat2;
use crate::path_tracing::{DirChange, Ray, RayInit, Ssp};
use num::complex::Complex64;
use num::Complex;
use pyo3::prelude::*;
use rayon::prelude::*;

use std::f64::consts::PI;

const TWO_PI: f64 = 2.0 * PI;
const I: Complex<f64> = Complex::new(0.0, 1.0);

/// Stores Beam propagation
#[derive(Clone)]
pub struct Beam {
    pub central_ray: Ray,
    pub p_vals: Vec<Complex<f64>>,
    pub q_vals: Vec<Complex<f64>>,
    pub c_vals: Vec<f64>,
    pub ang_step: f64,
}

/// Interface to python as Complex<f64> has no pyo3 type for conversion
#[pyclass]
#[derive(Debug)]
pub struct PyBeam {
    #[pyo3(get, set)]
    pub central_ray: Ray,
    #[pyo3(get, set)]
    pub p_re: Vec<f64>,
    #[pyo3(get, set)]
    pub q_re: Vec<f64>,
    #[pyo3(get, set)]
    pub p_im: Vec<f64>,
    #[pyo3(get, set)]
    pub q_im: Vec<f64>,
}

/// Stores locations to calculate pressure and calculated values corresponding to locations
pub struct PressureField {
    pub locations: Vec<[f64; 2]>,
    pub pressures: Vec<Complex<f64>>,
}

impl PyBeam {
    /// Initialise [`PyBeam`] from [`Beam`]
    pub fn from_beam(beam: &Beam) -> PyBeam {
        // bad and not good
        let q_re: Vec<f64> = beam.q_vals.iter().map(|&q| q.re).collect();
        let q_im: Vec<f64> = beam.q_vals.iter().map(|&q| q.im).collect();
        let p_re: Vec<f64> = beam.p_vals.iter().map(|&p| p.re).collect();
        let p_im: Vec<f64> = beam.p_vals.iter().map(|&p| p.im).collect();
        PyBeam {
            central_ray: beam.central_ray.clone(),
            p_re,
            q_re,
            p_im,
            q_im,
        }
    }
}

impl Beam {
    /// Initialise [`Beam`] struct from configs
    fn init_from_configs(init_source: &RayInit, prog_config: &ProgConfig) -> Beam {
        let mut bm: Beam = Beam {
            central_ray: Ray::init_from_cfgs(init_source, prog_config),
            p_vals: vec![Complex::ZERO; prog_config.max_it + 1],
            q_vals: vec![Complex::ZERO; prog_config.max_it + 1],
            c_vals: vec![0.0; prog_config.max_it + 1],
            ang_step: init_source.ang_step,
        };
        bm.q_vals[0] = Complex::new(0.0, 1.0 / init_source.init_sound_speed);
        bm.p_vals[0] = Complex::new(1.0, 0.0);
        bm.c_vals[0] = init_source.init_sound_speed;
        bm
    }

    /// Initialise [`Beam`] struct from configs and complete beam tracing process
    fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        env_config: &EnvConfig,
        ssp: &Ssp,
        pq_solver: &SolverMethod,
    ) -> Beam {
        let mut beam: Beam = Beam::init_from_configs(init_source, prog_config);
        let mut c_i: f64;
        let mut c_i1: f64;
        let mut c_i2: f64;
        let mut c_im1: f64;
        let mut g_i: f64;
        let ang: f64 = init_source.init_ang;
        let mut depth_dir: f64 = ang.sin().signum();
        // set initial sound speed values
        c_im1 = ssp.interp_sound_speed(
            beam.central_ray.depth_vals[0] - depth_dir * prog_config.depth_step,
        );
        c_i = ssp.interp_sound_speed(beam.central_ray.depth_vals[0]);
        c_i1 = ssp.interp_sound_speed(
            beam.central_ray.depth_vals[0] + depth_dir * prog_config.depth_step,
        );
        c_i2 = ssp.interp_sound_speed(
            beam.central_ray.depth_vals[0] + 2.0 * depth_dir * prog_config.depth_step,
        );

        for i in 0..env_config.isospaces.len() {
            if env_config.isospaces[i].body.contains_point(
                &beam.central_ray.range_vals[0],
                &beam.central_ray.depth_vals[0],
            ) {
                beam.central_ray
                    .iso_trace(&env_config.isospaces[i], prog_config.max_it, &ang);
                beam.update_pq_iso(&env_config.isospaces[i], 0);
                beam.central_ray.truncate_ray();
                beam.c_vals = vec![env_config.isospaces[i].sound_speed; beam.central_ray.ray_iter];
                return beam;
            }
        }

        // while loop to iterate updates of ray
        while (beam.central_ray.ray_iter < prog_config.max_it - init_source.init_iter)
            && (init_source
                .range_lims
                .contains(&beam.central_ray.range_vals[beam.central_ray.ray_iter]))
        {
            beam.c_vals[beam.central_ray.ray_iter] = c_i;
            // calculate local sound speed gradient
            g_i = (c_i1 - c_i) / prog_config.depth_step;
            let arc_step: f64 = (((beam.central_ray.ray_param * c_i1).asin()
                - (beam.central_ray.ray_param * c_i).asin())
                / (g_i * beam.central_ray.ray_param))
                .abs();
            // Update p-q equations
            beam.update_pq(
                &c_i,
                &c_i1,
                &c_im1,
                &c_i2,
                &arc_step,
                &prog_config.depth_step,
                pq_solver,
            );
            // iterate depth step. Match statement to update ssp variables as required
            match beam.central_ray.update_iteration(
                &c_i,
                &c_i1,
                &g_i,
                &depth_dir,
                &prog_config.depth_step,
            ) {
                DirChange::KeepDir => {
                    c_im1 = c_i;
                    c_i = c_i1;
                    c_i1 = c_i2;
                }
                DirChange::ChangeDir => {
                    std::mem::swap(&mut c_im1, &mut c_i1);
                    depth_dir = -depth_dir;
                }
            };
            // update c_i+2
            c_i2 = ssp.interp_sound_speed(
                beam.central_ray.depth_vals[beam.central_ray.ray_iter + 1]
                    + 2.0 * depth_dir * prog_config.depth_step,
            );

            // Check for intersections on all bodies in simulation
            if let Some(reflect_ans) = env_config.check_all_body_reflections(&beam.central_ray) {
                // Update intersection in step
                beam.central_ray.update_intersection(&reflect_ans, ssp);
                // recalculate depth direction
                depth_dir = reflect_ans.ang.sin().signum();
                // Interpolate for next sound speed profile value
                c_i = ssp.interp_sound_speed(reflect_ans.depth);
                c_im1 =
                    ssp.interp_sound_speed(reflect_ans.depth - depth_dir * prog_config.depth_step);
                c_i1 =
                    ssp.interp_sound_speed(reflect_ans.depth + depth_dir * prog_config.depth_step);
                c_i2 = ssp.interp_sound_speed(
                    reflect_ans.depth + 2.0 * depth_dir * prog_config.depth_step,
                );
            }
            // step iter value for calculation step taken
            beam.central_ray.ray_iter += 1;
        }

        // truncate vectors to remove any wasted space
        beam.truncate_beam();
        beam
    }

    /// Update p-q ODE values in [`IsoSpace`] propagation
    fn update_pq_iso(&mut self, isospace: &IsoSpace, start_id: usize) {
        // TODO: Update theory to match this
        for i in start_id..self.central_ray.ray_iter {
            let dist_step: f64 = ((self.central_ray.range_vals[i + 1]
                - self.central_ray.range_vals[i])
                .powi(2)
                + (self.central_ray.depth_vals[i + 1] - self.central_ray.depth_vals[i]).powi(2))
            .sqrt();
            self.p_vals[i + 1] = self.p_vals[i];
            self.q_vals[i + 1] = isospace.sound_speed * dist_step
                + Complex::I * isospace.sound_speed.powi(2)
                    / (self.central_ray.frequency * PI * self.ang_step.powi(2));
        }
    }

    /// Update p-q ODE values
    // TODO: Handle errors in a sane way here at some point  Ok(())
    #[allow(clippy::too_many_arguments)]
    fn update_pq(
        &mut self,
        c_i: &f64,
        c_i1: &f64,
        c_im1: &f64,
        c_i2: &f64,
        arc_step: &f64,
        depth_step: &f64,
        method: &SolverMethod,
    ) {
        match method {
            SolverMethod::RungeKutta4 => {
                self.update_pq_rk4(c_i, c_i1, c_im1, c_i2, arc_step, depth_step);
            }
            SolverMethod::Radau3IA => {
                self.update_pq_radau3_ia(c_i, c_i1, c_im1, c_i2, arc_step, depth_step);
            }
            SolverMethod::BackEuler => {
                self.update_pq_back_euler(c_i, c_i1, c_i2, arc_step, depth_step);
            }
            SolverMethod::Radau3IIA => {
                self.update_pq_radau3_iia(c_i, c_i1, c_im1, c_i2, arc_step, depth_step);
            }
        }
    }

    /// RK4 solver
    fn update_pq_rk4(
        &mut self,
        c_i: &f64,
        c_i1: &f64,
        c_im1: &f64,
        c_i2: &f64,
        arc_step: &f64,
        depth_step: &f64,
    ) {
        let c_i1_2: f64 = (c_i1 + c_i) / 2.0;
        // RK4 Method for now
        //calculate 2nd derivate of SSP wrt normal for steps j, j+1/2 and j+1
        let c_nn_j: f64 =
            -self.central_ray.ray_param * c_i * (c_i1 - 2.0 * c_i + c_im1) / depth_step.powi(2);
        let c_nn_j1_2: f64 =
            -self.central_ray.ray_param * c_i * (c_i2 - c_i1 - c_i + c_im1) / depth_step.powi(2);
        let c_nn_j1: f64 =
            -self.central_ray.ray_param * c_i1 * (c_i2 - 2.0 * c_i1 + c_i) / depth_step.powi(2);
        // calculate k values
        let k_n1: [Complex<f64>; 2] = [
            c_i * self.p_vals[self.central_ray.ray_iter],
            -c_nn_j * self.q_vals[self.central_ray.ray_iter] / c_i.powi(2),
        ];
        let k_n2: [Complex<f64>; 2] = [
            c_i1_2 * (self.p_vals[self.central_ray.ray_iter] + 0.5 * k_n1[0]),
            -c_nn_j1_2 * (self.q_vals[self.central_ray.ray_iter] + 0.5 * k_n1[1]) / c_i1_2.powi(2),
        ];
        let k_n3: [Complex<f64>; 2] = [
            c_i1_2 * (self.p_vals[self.central_ray.ray_iter] + 0.5 * k_n2[0]),
            -c_nn_j1_2 * (self.q_vals[self.central_ray.ray_iter] + 0.5 * k_n2[1]) / c_i1_2.powi(2),
        ];
        let k_n4: [Complex<f64>; 2] = [
            c_i1 * (self.p_vals[self.central_ray.ray_iter] + k_n3[0]),
            -c_nn_j1 * (self.q_vals[self.central_ray.ray_iter] + k_n3[1]) / c_i1.powi(2),
        ];
        // update p and q values
        self.q_vals[self.central_ray.ray_iter + 1] = self.q_vals[self.central_ray.ray_iter]
            + arc_step * (k_n1[0] + 2.0 * (k_n2[0] + k_n3[0]) + k_n4[0]) / 6.0;
        self.p_vals[self.central_ray.ray_iter + 1] = self.p_vals[self.central_ray.ray_iter]
            + arc_step * (k_n1[1] + 2.0 * (k_n2[1] + k_n3[1]) + k_n4[1]) / 6.0;
    }

    fn update_pq_radau3_iia(
        &mut self,
        c_i: &f64,
        c_i1: &f64,
        c_im1: &f64,
        c_i2: &f64,
        arc_step: &f64,
        depth_step: &f64,
    ) {
        let c_i1_3: f64 = (c_i1 + c_i + c_i) / 3.0;
        let c_nn_i: f64 =
            -self.central_ray.ray_param * c_i * (c_i1 - c_i - c_i + c_im1) / depth_step.powi(2);
        let c_nn_i1: f64 =
            -self.central_ray.ray_param * c_i1 * (c_i2 - c_i1 - c_i1 + c_i) / depth_step.powi(2);
        let c_nn_i1_3: f64 = (c_nn_i1 + c_nn_i + c_nn_i) / 3.0;
        let mat_ai1: Mat2<f64> = Mat2 {
            a: 0_f64,
            b: *c_i1,
            c: -c_nn_i1 / c_i1.powi(2),
            d: 0_f64,
        };
        let mat_ai1_3: Mat2<f64> = Mat2 {
            a: 0_f64,
            b: c_i1_3,
            c: -c_nn_i1_3 / c_i1_3.powi(2),
            d: 0_f64,
        };
        let coeff_mat: Mat2<f64> = (Mat2::I
            - (arc_step / 12_f64)
                * (5_f64 * mat_ai1_3 - 2_f64 * mat_ai1_3 * mat_ai1 + 3_f64 * mat_ai1))
            .inv()
            * (Mat2::I + (arc_step / 4_f64) * mat_ai1_3);
        let i: usize = self.central_ray.ray_iter;
        self.q_vals[i + 1] = coeff_mat.a * self.q_vals[i] + coeff_mat.b * self.p_vals[i];
        self.p_vals[i + 1] = coeff_mat.c * self.q_vals[i] + coeff_mat.d * self.p_vals[i];
    }

    /// Radau3IA Solver
    fn update_pq_radau3_ia(
        &mut self,
        c_i: &f64,
        c_i1: &f64,
        c_im1: &f64,
        c_i2: &f64,
        arc_step: &f64,
        depth_step: &f64,
    ) {
        let c_i2_3: f64 = (c_i1 + c_i1 + c_i) / 3.0;
        let c_nn_i: f64 =
            -self.central_ray.ray_param * c_i * (c_i1 - 2.0 * c_i + c_im1) / depth_step.powi(2);
        let c_nn_i1: f64 =
            -self.central_ray.ray_param * c_i1 * (c_i2 - 2.0 * c_i1 + c_i) / depth_step.powi(2);
        let c_nn_i2_3: f64 = (c_nn_i1 + c_nn_i1 + c_nn_i) / 3.0;
        let denom_i: f64 = 25_f64 * c_nn_i2_3 + c_i2_3;
        let cmat_i_11: f64 = 1_f64 + 15_f64 * c_nn_i2_3 / denom_i;
        let cmat_i_12: f64 =
            4_f64 * c_i.powi(2) / (arc_step * c_nn_i1) - 3_f64 * c_i2_3.powi(2) / denom_i;
        let cmat_i_21: f64 = 4_f64 / (arc_step * c_i) + 3_f64 * c_nn_i2_3 / (c_i2_3 * denom_i);
        let cmat_i_22: f64 = 1_f64 + 15_f64 * c_nn_i2_3 / denom_i;
        let dmat_i_11: f64 = -1_f64 - 15_f64 * c_nn_i2_3 / denom_i;
        let dmat_i_12: f64 = 3_f64 * c_i2_3.powi(2) / denom_i;
        let dmat_i_21: f64 = -3_f64 * c_nn_i2_3 / (c_i2_3 * denom_i);
        let dmat_i_22: f64 = -1_f64 - 15_f64 * c_nn_i2_3 / denom_i;
        let coeff: f64 = 4_f64 / (arc_step * (cmat_i_11 * cmat_i_22 - cmat_i_12 * cmat_i_21));
        let q_j: Complex<f64> = self.q_vals[self.central_ray.ray_iter];
        let p_j: Complex<f64> = self.p_vals[self.central_ray.ray_iter];
        let kn1_1: Complex<f64> = coeff
            * (q_j * (cmat_i_22 * dmat_i_11 - cmat_i_12 * dmat_i_21)
                + p_j * (cmat_i_22 * dmat_i_12 - cmat_i_12 * dmat_i_22));
        let kn1_2: Complex<f64> = coeff
            * (q_j * (-cmat_i_21 * dmat_i_11 + cmat_i_11 * dmat_i_21)
                + p_j * (-cmat_i_21 * dmat_i_12 + cmat_i_11 * dmat_i_22));
        let kn2_1: Complex<f64> = (-5_f64 * c_nn_i2_3 * (12_f64 * q_j / arc_step + 3_f64 * kn1_1)
            + c_i2_3.powi(2) * (12_f64 * p_j / arc_step + 3_f64 * kn1_2))
            / denom_i;
        let kn2_2: Complex<f64> = (-c_nn_i2_3 * (12_f64 * q_j / arc_step + 3_f64 * kn1_1) / c_i2_3
            - 5_f64 * c_i2_3.powi(2) * (12_f64 * p_j / arc_step + 3_f64 * kn1_2))
            / denom_i;
        self.q_vals[self.central_ray.ray_iter + 1] =
            q_j + arc_step * (kn1_1 + 3_f64 * kn2_1) / 4_f64;
        self.p_vals[self.central_ray.ray_iter + 1] =
            p_j + arc_step * (kn1_2 + 3_f64 * kn2_2) / 4_f64;
    }

    /// Backward Euler solver
    fn update_pq_back_euler(
        &mut self,
        c_i: &f64,
        c_i1: &f64,
        c_i2: &f64,
        arc_step: &f64,
        depth_step: &f64,
    ) {
        let c_nn_i1: f64 =
            -self.central_ray.ray_param * c_i1 * (c_i2 - 2.0 * c_i1 + c_i) / depth_step.powi(2);
        let coeff: f64 = c_i1 / (1.0 - arc_step.powi(2) * c_nn_i1);
        self.q_vals[self.central_ray.ray_iter + 1] = coeff
            * (self.q_vals[self.central_ray.ray_iter]
                + arc_step * c_i1 * self.p_vals[self.central_ray.ray_iter]);
        self.p_vals[self.central_ray.ray_iter + 1] = coeff
            * (self.p_vals[self.central_ray.ray_iter]
                - arc_step * c_nn_i1 * self.q_vals[self.central_ray.ray_iter])
            / c_i1.powi(2);
    }

    /// Remove unneeded empty cells
    pub fn truncate_beam(&mut self) {
        self.central_ray.truncate_ray();
        self.q_vals.truncate(self.central_ray.ray_iter + 1);
        self.p_vals.truncate(self.central_ray.ray_iter + 1);
    }

    fn calculate_pressure_predicates(
        &self,
        range_receiver: &f64,
        depth_receiver: &f64,
    ) -> Vec<PressurePredicate> {
        let mut predicates: Vec<PressurePredicate> = vec![];
        let mut range_step: f64;
        let mut depth_step: f64;
        let mut rng_intsct_perp: f64;
        let mut p: Complex<f64>;
        let mut q: Complex<f64>;
        let mut time: f64;
        let mut sound_speed: f64;
        let mut interp_frac: f64;
        let mut n: f64;

        for i in 0..self.central_ray.ray_iter {
            // TODO: move this function to math_util.rs as
            // interpolate_perpendicular_location
            range_step = self.central_ray.range_vals[i + 1] - self.central_ray.range_vals[i];
            depth_step = self.central_ray.depth_vals[i + 1] - self.central_ray.depth_vals[i];
            rng_intsct_perp = (range_step.powi(2) * range_receiver
                + depth_step.powi(2) * self.central_ray.range_vals[i]
                + range_step * depth_step * (depth_receiver - self.central_ray.depth_vals[i]))
                / (range_step.powi(2) + depth_step.powi(2));
            if ((rng_intsct_perp <= self.central_ray.range_vals[i + 1])
                && (rng_intsct_perp > self.central_ray.range_vals[i]))
                || ((rng_intsct_perp > self.central_ray.range_vals[i + 1])
                    && (rng_intsct_perp <= self.central_ray.range_vals[i]))
            {
                // TODO: MAke this an interp function????
                // FIXME: this wont work for vertical rays right now
                interp_frac = (rng_intsct_perp - self.central_ray.range_vals[i]) / range_step;
                n = ((range_receiver - self.central_ray.range_vals[i]) * depth_step
                    - (depth_receiver - self.central_ray.depth_vals[i]) * range_step)
                    / (range_step.powi(2) + depth_step.powi(2)).sqrt();
                p = self.p_vals[i] + interp_frac * (self.p_vals[i + 1] - self.p_vals[i]);
                q = self.q_vals[i] + interp_frac * (self.q_vals[i + 1] - self.q_vals[i]);
                time = self.central_ray.time_vals[i]
                    + interp_frac
                        * (self.central_ray.time_vals[i + 1] - self.central_ray.time_vals[i]);
                sound_speed = self.c_vals[i] + interp_frac * (self.c_vals[i + 1] - self.c_vals[i]);
                predicates.push(PressurePredicate {
                    p,
                    q,
                    time,
                    sound_speed,
                    n,
                });
            }
        }
        vec![]
    }

    pub fn calculate_pressure(&self, range: &f64, depth: &f64) -> Complex<f64> {
        let pre_pressures: Vec<PressurePredicate> =
            self.calculate_pressure_predicates(range, depth);
        let ang_freq: f64 = TWO_PI * self.central_ray.frequency;
        pre_pressures
            .iter()
            .map(|ent| ent.calculate_contribution(ang_freq))
            .sum::<Complex<f64>>()
    }
}

/// struct to hold required values for calculating pressure contribution of beam
struct PressurePredicate {
    p: Complex<f64>,
    q: Complex<f64>,
    time: f64,
    sound_speed: f64,
    n: f64,
}

impl PressurePredicate {
    /// Calculate ray coordinate contribution, disregarding spreading law
    pub fn calculate_contribution(&self, ang_freq: f64) -> Complex<f64> {
        (self.sound_speed / self.q).sqrt()
            * (I * ang_freq * (self.time + self.p * self.n.powi(2) / (self.q + self.q))).exp()
    }
}

pub fn evaluate_field(beam_config: BeamConfigRust, beams: Vec<Beam>) -> PressureField {
    // TODO: Check efficiency gains with par_iter() calls in different nest levels
    let pressures: Vec<Complex<f64>> = vec![I; beam_config.pressure_locs.len()];
    for range_depth in &beam_config.pressure_locs {
        beams
            .iter()
            .map(|bm| bm.calculate_pressure(&range_depth[0], &range_depth[1]))
            .sum::<Complex<f64>>();
    }
    PressureField {
        locations: beam_config.pressure_locs,
        pressures,
    }
}

/// Driving beam tracing function
pub fn trace_beams(cfg: BeamConfigRust) -> Vec<Beam> {
    let mut init_sources: Vec<RayInit> = vec![];
    let mut init_sound_speed: f64;
    for source in cfg.ray_config.sources {
        init_sound_speed = cfg
            .ray_config
            .env_config
            .ssp
            .interp_sound_speed(source.depth_pos);
        init_sources.append(
            &mut (0..source.n_rays)
                .map(|i| {
                    RayInit::from_source(&source, i, init_sound_speed, &cfg.ray_config.prog_config)
                })
                .collect(),
        );
    }
    init_sources
        .par_iter()
        .map(|init_source| {
            Beam::trace_from_init_source(
                init_source,
                &cfg.ray_config.prog_config,
                &cfg.ray_config.env_config,
                &cfg.ray_config.env_config.ssp,
                &cfg.pq_solver,
            )
        })
        .collect()
}
