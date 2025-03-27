use crate::interface::{Config, EnvConfig, ProgConfig};
use crate::path_tracing::{Ray, RayInit, Ssp, REFLECT_OFFSET};
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
    // pub p: Complex64,
    // #[pyo3(get, set)]
    // pub q: Complex64,
}

impl Beam {
    fn init_from_configs(init_source: &RayInit, prog_config: &ProgConfig) -> Beam {
        Beam {
            central_ray: Ray::init_from_cfgs(init_source, prog_config),
            // p: Complex64::new(1.0, 0.0),
            // q: Complex64::new(0.0, 1.0 / init_source.init_sound_speed),
        }
    }

    fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        env_config: &EnvConfig,
        ssp: &Ssp,
    ) -> Beam {
        let mut beam: Beam = Beam::init_from_configs(init_source, prog_config);
        let mut c_i: f64;
        let mut c_i1: f64;
        let mut g_i: f64;
        let ang: f64 = init_source.init_ang;
        let mut depth_dir: f64 = ang.sin().signum();
        // let mut range_dir: f64 = ang.cos().signum();
        // set initial sound speed values
        c_i = ssp.interp_sound_speed(beam.central_ray.depth_vals[0]);

        while (beam.central_ray.ray_iter < prog_config.max_it - init_source.init_iter)
            && (init_source
                .range_lims
                .contains(&beam.central_ray.range_vals[beam.central_ray.ray_iter]))
        {
            // set next SSP value
            c_i1 = ssp.interp_sound_speed(
                beam.central_ray.depth_vals[beam.central_ray.ray_iter]
                    + depth_dir * prog_config.depth_step,
            );
            // calculate local sound speed gradient
            g_i = (c_i1 - c_i) / prog_config.depth_step;
            // iterate depth step. This function will update value of c_i and depth_dir as required
            beam.central_ray.update_iteration(
                &mut c_i,
                &c_i1,
                &g_i,
                &mut depth_dir,
                &prog_config.depth_step,
            );
            // Check for intersections on all bodies in simulation
            if let Some(reflect_ans) = env_config.check_all_body_reflections(&beam.central_ray) {
                // Update intersection in step
                beam.central_ray.update_intersection(&reflect_ans, ssp);
                // recalculate depth direction
                depth_dir = reflect_ans.ang.sin().signum();
                // set next range and depth vals to minor offset from body to avoid ray getting
                // caught inside of [`Body`] struct
                beam.central_ray.range_vals[beam.central_ray.ray_iter + 1] =
                    reflect_ans.range + REFLECT_OFFSET * reflect_ans.ang.cos();
                beam.central_ray.depth_vals[beam.central_ray.ray_iter + 1] =
                    reflect_ans.depth + REFLECT_OFFSET * reflect_ans.ang.sin();
                // Interpolate for next sound speed profile value
                c_i = ssp.interp_sound_speed(reflect_ans.depth);
            }
            // step iter value for calculation step taken
            beam.central_ray.ray_iter += 1;
        }

        // save all rays to csv if config requires. NOTE:This is very slow
        if prog_config.save_to_csv {
            match beam.central_ray.write_to_csv(&prog_config.output_path) {
                Ok(_) => (),
                Err(_) => panic!(
                    "Ray {:} could not be written to csv",
                    beam.central_ray.ray_id
                ),
            }
        }

        // truncate vectors to remove any wasted space
        beam.central_ray
            .depth_vals
            .truncate(beam.central_ray.ray_iter + 1);
        beam.central_ray
            .range_vals
            .truncate(beam.central_ray.ray_iter + 1);
        beam
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
