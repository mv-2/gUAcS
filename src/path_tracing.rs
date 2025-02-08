use crate::interface::{Config, ProgConfig, SourceConfig};
use crate::math_util::deboor_alg;
use core::f64;
use csv::Writer;
use pyo3::prelude::*;
use rayon::prelude::*;
use serde::Deserialize;
use std::error::Error;
use std::ops::RangeInclusive;
use uuid::Uuid;

/// Stores data on body geometry
#[derive(Deserialize, Debug)]
pub struct Body {
    // polygonal definition only
    range_vals: Vec<f64>,
    depth_vals: Vec<f64>,
    // TODO: bounding box values to be added as later optimisation
}

/// Stores halfspace geometry and physical property data
#[derive(Deserialize, Debug)]
pub struct HalfSpace {
    body: Body,
    sound_speed: f64,
    density: f64,
}

/// Stores Ray propagation data
#[pyclass]
#[derive(Debug)]
pub struct Ray {
    #[pyo3(get, set)]
    pub range_vals: Vec<f64>,
    #[pyo3(get, set)]
    pub depth_vals: Vec<f64>,
    #[pyo3(get, set)]
    pub time_vals: Vec<f64>,
    #[pyo3(get, set)]
    pub ray_param: f64,
    #[pyo3(get, set)]
    pub ray_iter: usize,
    #[pyo3(get, set)]
    pub ray_id: String,
}

/// Stores data required to initialise rays
pub struct RayInit {
    range_pos: f64,
    depth_pos: f64,
    init_time: f64,
    init_ang: f64,
    init_iter: usize,
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
            init_iter: 0_usize,
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
        deboor_alg(depth, &self.ssp_knots, &self.ssp_coefs, self.ssp_degree)
    }
}

impl Ray {
    /// Trace ray using geometric theory
    pub fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        ssp: &Ssp,
        output_dir: &String,
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
            ray_iter: 0,
            ray_id: Uuid::new_v4().to_string(),
        };

        ray.range_vals[0] = init_source.range_pos;
        ray.depth_vals[0] = init_source.depth_pos;
        ray.time_vals[0] = init_source.init_time;

        c_i = ssp.interp_sound_speed(ray.depth_vals[0]);
        c_i1 = ssp.interp_sound_speed(ray.depth_vals[0] + depth_dir * prog_config.depth_step);

        while (ray.ray_iter < prog_config.max_it - init_source.init_iter)
            && (ray.range_vals[ray.ray_iter].abs() < prog_config.max_range)
        {
            g_i = (c_i1 - c_i) / prog_config.depth_step;
            if ray.ray_param * c_i1 < 1.0 {
                range_step = ((1.0 - (ray.ray_param * c_i).powi(2)).sqrt()
                    - (1.0 - (ray.ray_param * c_i1).powi(2)).sqrt())
                    / (ray.ray_param * g_i);
                time_step = ((c_i1 / c_i) * (1.0 + (1.0 - (ray.ray_param * c_i).powi(2)).sqrt())
                    / (1.0 + (1.0 - (ray.ray_param * c_i1).powi(2)).sqrt()))
                .ln()
                    / g_i.abs();
                ray.depth_vals[ray.ray_iter + 1] =
                    ray.depth_vals[ray.ray_iter] + depth_dir * prog_config.depth_step;
                c_i = c_i1;
            } else {
                ray.depth_vals[ray.ray_iter + 1] = ray.depth_vals[ray.ray_iter];
                depth_dir = -depth_dir;
                range_step = 2.0 * (1.0 - (ray.ray_param * c_i).powi(2)) / (ray.ray_param * g_i);
                time_step = 2.0
                    * ((1.0 + (1.0 - (ray.ray_param * c_i).powi(2)).sqrt())
                        / (ray.ray_param * c_i))
                        .ln()
                    / g_i.abs();
            }
            c_i1 = ssp.interp_sound_speed(
                ray.depth_vals[ray.ray_iter + 1] + depth_dir * prog_config.depth_step,
            );
            ray.range_vals[ray.ray_iter + 1] =
                ray.range_vals[ray.ray_iter] + range_dir * range_step;
            ray.time_vals[ray.ray_iter + 1] = ray.time_vals[ray.ray_iter] + time_step;
            ray.ray_iter += 1;
        }

        if prog_config.save_to_csv {
            match ray.write_to_csv(output_dir) {
                Ok(_) => (),
                Err(_) => panic!("Ray {:} could not be written to csv", ray.ray_id),
            }
        }
        ray.depth_vals.truncate(ray.ray_iter + 1);
        ray.range_vals.truncate(ray.ray_iter + 1);
        ray
    }

    pub fn write_to_csv(&self, output_dir: &String) -> Result<(), Box<dyn Error>> {
        let mut wtr = Writer::from_path(format!("{output_dir}/{:}.csv", self.ray_id))?;
        for i in 0..self.ray_iter {
            wtr.write_record(&[
                self.range_vals[i].to_string(),
                self.depth_vals[i].to_string(),
                self.time_vals[i].to_string(),
            ])?;
        }
        Ok(())
    }
}

impl Body {
    /// Checks for intersection between ['Ray'] and ['Body'] self struct over current ray
    /// interation step
    pub fn check_finite_intersection(&self, ray: &Ray) -> Option<(f64, f64, usize)> {
        let mut ray_dist_vals: Vec<Option<f64>> = vec![None; self.range_vals.len() - 1];
        let mut edge_dist: f64;
        let ray_range_step: f64 = ray.range_vals[ray.ray_iter + 1] - ray.range_vals[ray.ray_iter];
        let ray_depth_step: f64 = ray.depth_vals[ray.ray_iter + 1] - ray.depth_vals[ray.ray_iter];
        let mut edge_depth_step: f64;
        let mut edge_range_step: f64;
        let edge_range: RangeInclusive<f64> = RangeInclusive::new(0.0, 1.0);
        // iterate over each edge in polygon to find valid intersections
        for i in 0..(self.range_vals.len() - 1) {
            edge_depth_step = self.depth_vals[i + 1] - self.depth_vals[i];
            edge_range_step = self.range_vals[i + 1] - self.range_vals[i];
            // check for parallel edge and ray direction to save computation effort
            if edge_range_step / edge_depth_step != ray_range_step / ray_depth_step {
                edge_dist = ((self.range_vals[i] - ray.range_vals[ray.ray_iter]) * ray_depth_step
                    - (self.depth_vals[i] - ray.depth_vals[ray.ray_iter]) * ray_range_step)
                    / ((self.depth_vals[i + 1] - self.depth_vals[i]) * ray_range_step
                        - (self.range_vals[i + 1] - self.range_vals[i]) * ray_depth_step);
                // check if edge distance parameter is valid
                if edge_range.contains(&edge_dist) {
                    ray_dist_vals[i] = match ray_range_step >= 0.0 {
                        true => intersection_ray_dist_param(
                            edge_dist,
                            edge_range_step,
                            ray_range_step,
                            self.range_vals[i],
                            ray.range_vals[ray.ray_iter],
                        ),
                        false => intersection_ray_dist_param(
                            edge_dist,
                            edge_depth_step,
                            ray_depth_step,
                            self.depth_vals[i],
                            ray.depth_vals[ray.ray_iter],
                        ),
                    }
                }
            }
        }
        // Check if there are any valid intersections calculated
        if !ray_dist_vals
            .iter()
            .filter(|&x| !x.is_none())
            .collect::<Vec<&Option<f64>>>()
            .is_empty()
        {
            let ray_dists_numeric: Vec<f64> = ray_dist_vals
                .iter()
                .map(|&x| x.unwrap_or(f64::INFINITY))
                .collect();
            // Using direct unwrap here as None values should not be present in this and a minimum
            // value should always exist
            let edge_id: usize = ray_dists_numeric
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(id, _)| id)
                .expect(
                    "There is no minimum distance when calculating intersection for {ray.ray_id}",
                );
            Some((
                ray.range_vals[ray.ray_iter] + ray_dists_numeric[edge_id] * ray_range_step,
                ray.depth_vals[ray.ray_iter] + ray_dists_numeric[edge_id] * ray_depth_step,
                edge_id,
            ))
        } else {
            None
        }
    }

    pub fn calculate_reflection(&self, ray: &Ray) -> Option<f64> {
        // let intersection_ans: (f64, f64, usize) = match self.check_finite_intersection(ray) {
        //     Some(ans) => ans,
        //     None => return None,
        // };
        let intersection_ans: (f64, f64, usize) = self.check_finite_intersection(ray)?;
        //TODO: Add time interpolation
        // ray.range_vals[ray.ray_iter + 1] = intersection_ans.0;
        // ray.depth_vals[ray.ray_iter + 1] = intersection_ans.1;
        let ray_ang: f64 = (ray.depth_vals[ray.ray_iter + 1] - ray.depth_vals[ray.ray_iter])
            .atan2(ray.range_vals[ray.ray_iter + 1] - ray.range_vals[ray.ray_iter]);
        let side_ang: f64 = self.get_edge_ang(intersection_ans.2);
        Some(ray_ang + 2.0 * side_ang)
    }

    // calculates angle of edge # edge_id
    //TODO: check wrapping
    fn get_edge_ang(&self, edge_id: usize) -> f64 {
        (self.depth_vals[edge_id + 1] - self.depth_vals[edge_id])
            .atan2(self.range_vals[edge_id + 1] - self.range_vals[edge_id])
    }
}

/// function is range/depth agnostic as detailed in theory document
fn intersection_ray_dist_param(
    edge_dist_param: f64,
    edge_step: f64,
    ray_step: f64,
    edge_point_0: f64,
    ray_point_0: f64,
) -> Option<f64> {
    let dist: f64 = ((edge_point_0 - ray_point_0) + edge_dist_param * edge_step) / ray_step;
    if dist.is_sign_positive() {
        Some(dist)
    } else {
        None
    }
}

/// Driver function to trace all rays from [`Config`] struct
pub fn trace_from_config(cfg: Config, output_dir: String) -> Vec<Ray> {
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
            Ray::trace_from_init_source(
                ray_source,
                &cfg.prog_config,
                &cfg.env_config.ssp,
                &output_dir,
            )
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn ray_intersection_tests() {
        // test 1 expect success
        let test_body: Body = Body {
            range_vals: vec![-1.0, 1.0, 1.0, -1.0, -1.0],
            depth_vals: vec![-1.0, -1.0, 1.0, 1.0, -1.0],
        };

        let test_ray: Ray = Ray {
            range_vals: vec![-10.0, 10.0],
            depth_vals: vec![-1.0, 1.0],
            ray_iter: 0,
            ray_id: Uuid::new_v4().to_string(),
            ray_param: 1.0,
            time_vals: vec![0.0, 1.0],
        };
        let ans: (f64, f64, usize) = test_body.check_finite_intersection(&test_ray).unwrap();
        assert!(ans.2 == 3);
        assert!((ans.1 + 1.0 - 9.0 / 10.0).abs() < 1e-8);
        assert!((ans.0 + 1.0).abs() < 1e-8);

        // test 2 parrallel line failure
        let test_body: Body = Body {
            range_vals: vec![0.0, 100.0],
            depth_vals: vec![0.0, 0.0],
        };

        let test_ray: Ray = Ray {
            range_vals: vec![-10.0, 10.0],
            depth_vals: vec![-1.0, -1.0],
            ray_iter: 0,
            ray_id: Uuid::new_v4().to_string(),
            ray_param: 1.0,
            time_vals: vec![0.0, 1.0],
        };
        let ans: Option<(f64, f64, usize)> = test_body.check_finite_intersection(&test_ray);
        assert!(ans.is_none());
    }

    #[test]
    fn reflection_test() {
        // test 1 expect success
        let test_body: Body = Body {
            range_vals: vec![-1.0, 1.0, 1.0, -1.0, -1.0],
            depth_vals: vec![-1.0, -1.0, 1.0, 1.0, -1.0],
        };

        let test_ray: Ray = Ray {
            range_vals: vec![-10.0, 0.0],
            depth_vals: vec![-1.0, 0.0],
            ray_iter: 0,
            ray_id: Uuid::new_v4().to_string(),
            ray_param: 1.0,
            time_vals: vec![0.0, 1.0],
        };
        // TODO: VERIFY THIS
        let ang_out: f64 = test_body.calculate_reflection(&test_ray).unwrap();
        println!("{ang_out}");
    }
}
