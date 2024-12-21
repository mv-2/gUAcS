use crate::interface::{Config, ProgConfig, SourceConfig};
use crate::math_util::deboor_alg;
use rayon::prelude::*;
use serde::Deserialize;

/// Stores data on body geometry
#[derive(Deserialize, Debug)]
pub struct Body {
    // polygonal definition only
    range_vals: Vec<f64>,
    depth_vals: Vec<f64>,
    // TODO bounding box values to be added as later optimisation
}

/// Stores halfspace geometry and physical property data
#[derive(Deserialize, Debug)]
pub struct HalfSpace {
    body: Body,
    sound_speed: f64,
    density: f64,
}

/// Stores Ray propagation data
#[derive(Debug)]
pub struct Ray {
    pub range_vals: Vec<f64>,
    pub depth_vals: Vec<f64>,
    pub time_vals: Vec<f64>,
    pub ray_param: f64,
}

/// Stores data required to initialise rays
pub struct RayInit {
    range_pos: f64,
    depth_pos: f64,
    init_time: f64,
    init_ang: f64,
}

impl RayInit {
    /// Create [`RayInit`] struct from [`SourceConfig`] struct
    fn from_source(source: &SourceConfig, ray_id: usize) -> Self {
        let init_ang: f64 = source.ray_fan_limits[0]
            + (ray_id as f64) * (source.ray_fan_limits[1] - source.ray_fan_limits[0])
                / (source.n_rays as f64 - 1_f64);

        RayInit {
            range_pos: source.range_pos,
            depth_pos: source.depth_pos,
            init_time: 0.0,
            init_ang,
        }
    }
}

/// Stores sound speed profile data
#[derive(Deserialize, Debug)]
pub struct Ssp {
    ssp_knots: Vec<f64>,
    ssp_coefs: Vec<f64>,
    ssp_degree: usize,
}

impl Ssp {
    /// Calculate sound speed value at given depth
    pub fn interp_sound_speed(&self, depth: f64) -> f64 {
        deboor_alg(depth, &self.ssp_knots, &self.ssp_coefs, &self.ssp_degree)
    }
}

impl Ray {
    /// Trace ray using geometric theory
    pub fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        ssp: &Ssp,
    ) -> Self {
        let mut range_step: f64;
        let mut time_step: f64;
        let mut c_i: f64;
        let mut c_i1: f64;
        let mut g_i: f64;
        let ang: f64 = init_source.init_ang;
        let mut depth_dir: f64 = ang.sin().signum();
        let range_dir: f64 = ang.cos().signum();

        let mut ray: Ray = Ray {
            range_vals: vec![0.0; prog_config.max_it + 1],
            depth_vals: vec![0.0; prog_config.max_it + 1],
            time_vals: vec![0.0; prog_config.max_it + 1],
            ray_param: ang.cos() / ssp.interp_sound_speed(init_source.depth_pos),
        };

        ray.range_vals[0] = init_source.range_pos;
        ray.depth_vals[0] = init_source.depth_pos;
        ray.time_vals[0] = init_source.init_time;

        for i in 0..prog_config.max_it {
            c_i = ssp.interp_sound_speed(ray.depth_vals[i]);
            c_i1 = ssp.interp_sound_speed(ray.depth_vals[i] + depth_dir * prog_config.depth_step);
            g_i = (c_i1 - c_i) / 2.0;
            if (ray.ray_param * c_i1).powi(2) < 1.0 {
                range_step = ((1.0 - (ray.ray_param * c_i).powi(2)).sqrt()
                    - (1.0 - (ray.ray_param * c_i1).powi(2)).sqrt())
                    / (ray.ray_param * g_i);
                time_step = ((c_i1 / c_i) * (1.0 + (1.0 - (ray.ray_param * c_i).powi(2)).sqrt())
                    / (1.0 + (1.0 - (ray.ray_param * c_i1).powi(2)).sqrt()))
                .ln()
                    / g_i.abs();
                ray.depth_vals[i + 1] = ray.depth_vals[i] + depth_dir * prog_config.depth_step;
            } else {
                depth_dir = -depth_dir;
                range_step = 2.0 * (1.0 - (ray.ray_param * c_i).powi(2)) / (ray.ray_param * g_i);
                time_step = 2.0
                    * ((1.0 + (1.0 - (ray.ray_param * c_i).powi(2)).sqrt())
                        / (ray.ray_param * c_i))
                        .ln()
                    / g_i.abs();
            }
            ray.range_vals[i + 1] = ray.range_vals[i] + range_dir * range_step;
            ray.time_vals[i + 1] = ray.time_vals[i] + time_step;
        }
        ray
    }
}

pub fn trace_from_config(cfg: Config) -> Vec<Ray> {
    let mut init_sources: Vec<RayInit> = vec![];
    for source in cfg.sources {
        init_sources.append(
            &mut (0..source.n_rays)
                .map(|i| RayInit::from_source(&source, i))
                .collect(),
        );
    }
    init_sources
        .par_iter()
        .map(|ray_source| {
            Ray::trace_from_init_source(ray_source, &cfg.prog_config, &cfg.env_config.ssp)
        })
        .collect()
}
