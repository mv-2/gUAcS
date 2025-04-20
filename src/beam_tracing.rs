use crate::interface::{Config, EnvConfig, ProgConfig};
use crate::path_tracing::{DirChange, Ray, RayInit, Ssp};
use num::complex::Complex64;
use pyo3::prelude::*;
use rayon::prelude::*;

// Stores Beam propagation
#[pyclass]
#[derive(Debug)]
pub struct Beam {
    #[pyo3(get, set)]
    pub central_ray: Ray,
    // #[pyo3(get, set)]
    pub p_vals: Vec<Complex64>,
    // #[pyo3(get, set)]
    pub q_vals: Vec<Complex64>,
    // #[pyo3(get, set)]
    pub eps_vals: Vec<Complex64>,
}

impl Beam {
    /// Initialise [`Beam`] struct from configs
    fn init_from_configs(init_source: &RayInit, prog_config: &ProgConfig) -> Beam {
        let mut bm: Beam = Beam {
            central_ray: Ray::init_from_cfgs(init_source, prog_config),
            p_vals: vec![Complex64::ZERO; prog_config.max_it + 1],
            q_vals: vec![Complex64::ZERO; prog_config.max_it + 1],
            eps_vals: vec![Complex64::ZERO; prog_config.max_it + 1],
        };
        bm.q_vals[0] = Complex64::new(0.0, 1.0 / init_source.init_sound_speed);
        bm.q_vals[0] = Complex64::new(1.0, 0.0);
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

        while (beam.central_ray.ray_iter < prog_config.max_it - init_source.init_iter)
            && (init_source
                .range_lims
                .contains(&beam.central_ray.range_vals[beam.central_ray.ray_iter]))
        {
            // calculate local sound speed gradient
            g_i = (c_i1 - c_i) / prog_config.depth_step;
            // Update p-q equations
            beam.update_pq(&c_i, &c_i1, &c_im1, &c_i2, &g_i, &prog_config.depth_step);
            // iterate depth step. Match statement to update ssp variables as required
            match beam.central_ray.update_iteration(
                &c_i,
                &c_i1,
                &g_i,
                &mut depth_dir,
                &prog_config.depth_step,
            ) {
                DirChange::KeepDir => {
                    c_im1 = c_i;
                    c_i = c_i1;
                    c_i1 = c_i2;
                }
                DirChange::ChangeDir => {
                    std::mem::swap(&mut c_im1, &mut c_i1);
                }
            };
            c_i2 = ssp.interp_sound_speed(
                beam.central_ray.depth_vals[beam.central_ray.ray_iter]
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
    fn update_pq(
        &mut self,
        c_i: &f64,
        c_i1: &f64,
        c_im1: &f64,
        c_i2: &f64,
        g_i: &f64,
        depth_step: &f64,
    ) {
        let c_i1_2: f64 = (c_i1 + c_i) / 2.0;
        // Calculate arc length in last depth step
        let arc_step: f64 = ((self.central_ray.ray_param * c_i1 + 1.0)
            * (self.central_ray.ray_param * c_i - 1.0)
            / ((self.central_ray.ray_param * c_i1 - 1.0)
                * (self.central_ray.ray_param * c_i + 1.0)))
            .abs()
            .ln()
            / (2.0 * g_i * self.central_ray.ray_param);
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

    /// Remove unneeded empty cells
    pub fn truncate_beam(&mut self) {
        self.central_ray.truncate_ray();
        self.eps_vals.truncate(self.central_ray.ray_iter + 1);
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
