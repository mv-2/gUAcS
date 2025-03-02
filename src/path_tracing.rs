use crate::interface::{Config, EnvConfig, ProgConfig, SourceConfig};
use crate::math_util::deboor_alg;
use core::f64;
use csv::Writer;
use pyo3::prelude::*;
use rayon::prelude::*;
use serde::Deserialize;
use std::error::Error;
use std::ops::RangeInclusive;
use uuid::Uuid;

// NUCLEAR OPTION FOR REFLECTION CALCULATION
const REFECT_OFFSET: f64 = 0.1;

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

struct IntersectLoc {
    range: f64,
    depth: f64,
    edge_id: usize,
}

struct ReflectResult {
    range: f64,
    depth: f64,
    ang: f64,
}

impl Ssp {
    /// Calculate sound speed value at given depth
    pub fn interp_sound_speed(&self, depth: f64) -> f64 {
        deboor_alg(depth, &self.ssp_knots, &self.ssp_coefs, self.ssp_degree)
    }
}

impl EnvConfig {
    fn check_all_body_reflections(&self, ray: &Ray) -> Option<ReflectResult> {
        // let mut temp_ans: Option<(f64, f64, usize)>;
        // I would love to directly iterate body values but unfortunately borrow checker is whiny
        for i in 0..self.bodies.len() {
            match self.bodies[i].calculate_reflection(ray) {
                Some(ans) => return Some(ans),
                None => 1,
            };
            // temp_ans = self.bodies[i].check_finite_intersection(ray);
            // if temp_ans.is_some() {
            // return temp_ans;
            // }
        }
        None
    }
}

impl Ray {
    /// Trace ray using geometric theory
    pub fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        env_config: &EnvConfig,
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
        let mut range_dir: f64 = ang.cos().signum();

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
                range_step = (((1.0 - (ray.ray_param * c_i).powi(2)).sqrt()
                    - (1.0 - (ray.ray_param * c_i1).powi(2)).sqrt())
                    / (ray.ray_param * g_i));
                time_step = ((c_i1 / c_i) * (1.0 + (1.0 - (ray.ray_param * c_i).powi(2)).sqrt())
                    / (1.0 + (1.0 - (ray.ray_param * c_i1).powi(2)).sqrt()))
                .ln()
                    / g_i;
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
            // Update range
            ray.range_vals[ray.ray_iter + 1] =
                ray.range_vals[ray.ray_iter] + range_dir * range_step;
            ray.time_vals[ray.ray_iter + 1] = ray.time_vals[ray.ray_iter] + time_step;
            // Check if any intersections occur and if any valid solutions exist then ray
            // parameters will be reassigned
            if let Some(ans) = env_config.check_all_body_reflections(&ray) {
                ray.range_vals[ray.ray_iter + 1] = ans.range;
                ray.depth_vals[ray.ray_iter + 1] = ans.depth;
                ray.ray_param = ans.ang.cos() / ssp.interp_sound_speed(ans.depth); // TODO: Check if this makes any major
                                                                                   // difference of if extra interpolation can be afforded
                depth_dir = ans.ang.sin().signum();
                range_dir = ans.ang.cos().signum();
                ray.ray_iter += 1;
                ray.range_vals[ray.ray_iter + 1] = ans.range + REFECT_OFFSET * ans.ang.cos();
                ray.depth_vals[ray.ray_iter + 1] = ans.depth + REFECT_OFFSET * ans.ang.sin();
                c_i = ssp.interp_sound_speed(ans.depth);
            }
            // assign next ssp value
            c_i1 = ssp.interp_sound_speed(
                ray.depth_vals[ray.ray_iter + 1] + depth_dir * prog_config.depth_step,
            );
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
    ///
    /// TODO: Come back to this and see if logic can be compressed because this is too long for its
    /// functionality at this point
    #[allow(clippy::needless_range_loop)]
    fn check_finite_intersection(&self, ray: &Ray) -> Option<IntersectLoc> {
        let mut edge_dist: f64;
        let mut edge_depth_step: f64;
        let mut edge_range_step: f64;
        let mut ray_dist_vals: Vec<Option<f64>> = vec![None; self.range_vals.len() - 1];
        let ray_range_step: f64 = ray.range_vals[ray.ray_iter + 1] - ray.range_vals[ray.ray_iter];
        let ray_depth_step: f64 = ray.depth_vals[ray.ray_iter + 1] - ray.depth_vals[ray.ray_iter];
        // Define valid ranges for solution parameters as defined in theory document
        let edge_param_range: RangeInclusive<f64> = RangeInclusive::new(0.0, 1.0);
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
                // check if edge distance parameter is in valid range
                if edge_param_range.contains(&edge_dist) {
                    // ensure that solution is defined without zero denominator
                    ray_dist_vals[i] = match ray_range_step != 0.0 {
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
        // convert options to numeric results
        let ray_dists_numeric: Vec<f64> = ray_dist_vals
            .iter()
            .map(|&x| x.unwrap_or(f64::INFINITY))
            .collect();
        // calculate minimum distance
        let edge_id: Option<usize> = ray_dists_numeric
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(id, _)| id);
        // Check if valid minimiser id exists
        match edge_id {
            // if valid minimiser exists then check if hit distance is valid otherwise return None
            Some(id) => match edge_param_range.contains(&ray_dists_numeric[id]) {
                true => Some(IntersectLoc {
                    range: ray.range_vals[ray.ray_iter] + ray_dists_numeric[id] * ray_range_step,
                    depth: ray.depth_vals[ray.ray_iter] + ray_dists_numeric[id] * ray_depth_step,
                    edge_id: id,
                }),
                false => None,
            },
            None => None,
        }
    }

    fn calculate_reflection(&self, ray: &Ray) -> Option<ReflectResult> {
        let intersection_ans: IntersectLoc = self.check_finite_intersection(ray)?;
        //TODO: Add time interpolation
        let ray_ang: f64 = (ray.depth_vals[ray.ray_iter + 1] - ray.depth_vals[ray.ray_iter])
            .atan2(ray.range_vals[ray.ray_iter + 1] - ray.range_vals[ray.ray_iter]);
        let side_ang: f64 = self.get_edge_ang(intersection_ans.edge_id);
        let new_ang: f64 = 2.0 * side_ang - ray_ang;
        Some(ReflectResult {
            range: intersection_ans.range,
            depth: intersection_ans.depth,
            ang: new_ang,
        })
    }

    // calculates angle of edge # edge_id
    //TODO: check wrapping
    // Change this to eval on compile then use this function as a lookup for increased
    // performance??
    fn get_edge_ang(&self, edge_id: usize) -> f64 {
        let mut ang: f64 = (self.depth_vals[edge_id + 1] - self.depth_vals[edge_id])
            .atan2(self.range_vals[edge_id + 1] - self.range_vals[edge_id])
            % f64::consts::PI;
        ang = (ang + f64::consts::PI) % f64::consts::PI;
        match ang > f64::consts::PI / 2.0 {
            true => ang - f64::consts::PI,
            false => ang,
        }
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
                &cfg.env_config,
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
        let ans: IntersectLoc = test_body.check_finite_intersection(&test_ray).unwrap();
        assert!(ans.edge_id == 3);
        assert!((ans.range + 1.0 - 9.0 / 10.0).abs() < 1e-8);
        assert!((ans.depth + 1.0).abs() < 1e-8);

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
        let ans: Option<IntersectLoc> = test_body.check_finite_intersection(&test_ray);
        assert!(ans.is_none());

        // test 3 failure by angle
        let test_body: Body = Body {
            range_vals: vec![-1.0, 1.0, 1.0, -1.0, -1.0],
            depth_vals: vec![-1.0, -1.0, 1.0, 1.0, -1.0],
        };

        let test_ray: Ray = Ray {
            range_vals: vec![-10.0, 1.0],
            depth_vals: vec![-1.0, -120.0],
            ray_iter: 0,
            ray_id: Uuid::new_v4().to_string(),
            ray_param: 1.0,
            time_vals: vec![0.0, 1.0],
        };
        let ans: Option<IntersectLoc> = test_body.check_finite_intersection(&test_ray);
        assert!(ans.is_none());

        // test 4 failure by short
        let test_body: Body = Body {
            range_vals: vec![-1.0, 1.0, 1.0, -1.0, -1.0],
            depth_vals: vec![-1.0, -1.0, 1.0, 1.0, -1.0],
        };

        let test_ray: Ray = Ray {
            range_vals: vec![-10.0, -5.0],
            depth_vals: vec![-1.0, -0.5],
            ray_iter: 0,
            ray_id: Uuid::new_v4().to_string(),
            ray_param: 1.0,
            time_vals: vec![0.0, 1.0],
        };
        let ans: Option<IntersectLoc> = test_body.check_finite_intersection(&test_ray);
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
        let ang_out: ReflectResult = test_body.calculate_reflection(&test_ray).unwrap();
        println!("{:}", ang_out.ang);
    }
}
