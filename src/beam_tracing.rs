use crate::interface::{Config, EnvConfig, ProgConfig};
use crate::math_util::Mat2;
use crate::path_tracing::{DirChange, Ray, RayInit, Ssp};
use num::complex::Complex64;
use pyo3::prelude::*;
use rayon::prelude::*;

/// Stores Beam propagation
#[derive(Clone)]
pub struct Beam {
    pub central_ray: Ray,
    pub p_vals: Vec<Complex64>,
    pub q_vals: Vec<Complex64>,
}

/// Interface to python as Complex64 has no pyo3 type for conversion
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

// Enum to set solver method for p-q equations
enum SolverMethod {
    RungeKutta4,
    Radau3IA,
    BackEuler,
    Radau3IIA,
}
impl SolverMethod {
    pub fn from_string(str: &str) -> Option<SolverMethod> {
        match str {
            "BackwardEuler" => Some(SolverMethod::BackEuler),
            "RungeKutta4" => Some(SolverMethod::RungeKutta4),
            "Radau3IA" => Some(SolverMethod::Radau3IA),
            "Radau3IIA" => Some(SolverMethod::Radau3IIA),
            _ => None,
        }
    }
}

impl Beam {
    /// Initialise [`Beam`] struct from configs
    fn init_from_configs(init_source: &RayInit, prog_config: &ProgConfig) -> Beam {
        let mut bm: Beam = Beam {
            central_ray: Ray::init_from_cfgs(init_source, prog_config),
            p_vals: vec![Complex64::ZERO; prog_config.max_it + 1],
            q_vals: vec![Complex64::ZERO; prog_config.max_it + 1],
        };
        bm.q_vals[0] = Complex64::new(0.0, 1.0 / init_source.init_sound_speed);
        bm.p_vals[0] = Complex64::new(1.0, 0.0);
        bm
    }

    /// Initialise [`Beam`] struct from configs and complete beam tracing process
    fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        env_config: &EnvConfig,
        ssp: &Ssp,
    ) -> Beam {
        let mut beam: Beam = Beam::init_from_configs(init_source, prog_config);
        let method: SolverMethod = SolverMethod::from_string(&prog_config.pq_solver[..])
            .expect("Requested p-q solver method is not defined");
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

        // while loop to iterate updates of ray
        while (beam.central_ray.ray_iter < prog_config.max_it - init_source.init_iter)
            && (init_source
                .range_lims
                .contains(&beam.central_ray.range_vals[beam.central_ray.ray_iter]))
        {
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
                &method,
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
                    println!("{:}", beam.central_ray.ray_iter);
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
        let k_n1: [Complex64; 2] = [
            c_i * self.p_vals[self.central_ray.ray_iter],
            -c_nn_j * self.q_vals[self.central_ray.ray_iter] / c_i.powi(2),
        ];
        let k_n2: [Complex64; 2] = [
            c_i1_2 * (self.p_vals[self.central_ray.ray_iter] + 0.5 * k_n1[0]),
            -c_nn_j1_2 * (self.q_vals[self.central_ray.ray_iter] + 0.5 * k_n1[1]) / c_i1_2.powi(2),
        ];
        let k_n3: [Complex64; 2] = [
            c_i1_2 * (self.p_vals[self.central_ray.ray_iter] + 0.5 * k_n2[0]),
            -c_nn_j1_2 * (self.q_vals[self.central_ray.ray_iter] + 0.5 * k_n2[1]) / c_i1_2.powi(2),
        ];
        let k_n4: [Complex64; 2] = [
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
        let c_i2_3 = (c_i1 + c_i1 + c_i) / 3.0;
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
        let q_j: Complex64 = self.q_vals[self.central_ray.ray_iter];
        let p_j: Complex64 = self.p_vals[self.central_ray.ray_iter];
        let kn1_1: Complex64 = coeff
            * (q_j * (cmat_i_22 * dmat_i_11 - cmat_i_12 * dmat_i_21)
                + p_j * (cmat_i_22 * dmat_i_12 - cmat_i_12 * dmat_i_22));
        let kn1_2: Complex64 = coeff
            * (q_j * (-cmat_i_21 * dmat_i_11 + cmat_i_11 * dmat_i_21)
                + p_j * (-cmat_i_21 * dmat_i_12 + cmat_i_11 * dmat_i_22));
        let kn2_1: Complex64 = (-5_f64 * c_nn_i2_3 * (12_f64 * q_j / arc_step + 3_f64 * kn1_1)
            + c_i2_3.powi(2) * (12_f64 * p_j / arc_step + 3_f64 * kn1_2))
            / denom_i;
        let kn2_2: Complex64 = (-c_nn_i2_3 * (12_f64 * q_j / arc_step + 3_f64 * kn1_1) / c_i2_3
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
}

/// Driving beam tracing function
pub fn trace_beams(cfg: Config) -> Vec<Beam> {
    let mut init_sources: Vec<RayInit> = vec![];
    let mut init_sound_speed: f64;
    for source in cfg.sources {
        init_sound_speed = cfg.env_config.ssp.interp_sound_speed(source.depth_pos);
        init_sources.append(
            &mut (0..source.n_rays)
                .map(|i| RayInit::from_source(&source, i, init_sound_speed, &cfg.prog_config))
                .collect(),
        );
    }
    init_sources
        .par_iter()
        .map(|init_source| {
            Beam::trace_from_init_source(
                init_source,
                &cfg.prog_config,
                &cfg.env_config,
                &cfg.env_config.ssp,
            )
        })
        .collect()
}
