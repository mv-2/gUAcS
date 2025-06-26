use crate::interface::{Config, EnvConfig, IsoSpace, ProgConfig, SourceConfig};
use crate::math_util::deboor_alg;
use core::f64;
use csv::Writer;
use pyo3::prelude::*;
use rayon::prelude::*;
use serde::Deserialize;
use std::error::Error;
use std::ops::RangeInclusive;
use uuid::Uuid;

//TODO: NUCLEAR OPTION FOR REFLECTION CALCULATION
pub const REFLECT_OFFSET: f64 = 0.1;

/// Stores data on body geometry
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct Body {
    // polygonal definition only
    #[pyo3(get, set)]
    pub range_vals: Vec<f64>,
    #[pyo3(get, set)]
    pub depth_vals: Vec<f64>,
    // TODO: bounding box values to be added as later optimisation
}

///// Stores halfspace geometry and physical property data
//#[derive(Deserialize, Debug, Clone)]
//#[pyclass]
//pub struct HalfSpace {
//    #[pyo3(get, set)]
//    pub body: Body,
//    #[pyo3(get, set)]
//    pub sound_speed: f64,
//    #[pyo3(get, set)]
//    pub density: f64,
//}

/// Stores Ray propagation data
#[pyclass]
#[derive(Debug, Clone)]
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
    pub range_pos: f64,
    pub depth_pos: f64,
    pub init_time: f64,
    pub init_ang: f64,
    pub init_iter: usize,
    pub init_sound_speed: f64,
    pub range_lims: RangeInclusive<f64>,
    pub frequency: f64,
}

impl RayInit {
    /// Create [`RayInit`] struct from [`SourceConfig`] struct
    pub fn from_source(
        source: &SourceConfig,
        ray_id: usize,
        init_sound_speed: f64,
        prog_config: &ProgConfig,
    ) -> Self {
        let init_ang: f64 = match source.n_rays > 1 {
            true => {
                source.ray_fan_limits[0]
                    + (ray_id as f64) * (source.ray_fan_limits[1] - source.ray_fan_limits[0])
                        / (source.n_rays as f64 - 1_f64)
            }
            false => source.ray_fan_limits[0],
        };

        RayInit {
            range_pos: source.range_pos,
            depth_pos: source.depth_pos,
            init_time: 0.0,
            init_ang,
            init_iter: 0_usize,
            init_sound_speed,
            range_lims: RangeInclusive::new(prog_config.min_range, prog_config.max_range), // TODO:
            // Initialise this elsewhere
            frequency: source.frequency,
        }
    }
}

/// Stores sound speed profile data
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct Ssp {
    #[pyo3(get, set)]
    pub ssp_knots: Vec<f64>,
    #[pyo3(get, set)]
    pub ssp_coefs: Vec<f64>,
    #[pyo3(get, set)]
    pub ssp_degree: usize,
}

/// Stores result of [`check_finite_intersection()`]
struct IntersectLoc {
    range: f64,
    depth: f64,
    time: f64,
    edge_id: usize,
}

/// Stores result of reflection calculation in [`calculate_reflection()`]
pub struct ReflectResult {
    pub range: f64,
    pub depth: f64,
    pub time: f64,
    pub ang: f64,
}

impl Ssp {
    /// Calculate sound speed value at given depth
    pub fn interp_sound_speed(&self, depth: f64) -> f64 {
        deboor_alg(depth, &self.ssp_knots, &self.ssp_coefs, self.ssp_degree)
    }
}

impl EnvConfig {
    pub fn check_all_body_reflections(&self, ray: &Ray) -> Option<ReflectResult> {
        // let mut temp_ans: Option<(f64, f64, usize)>;
        // I would love to directly iterate body values but unfortunately borrow checker is whiny
        for i in 0..self.bodies.len() {
            match self.bodies[i].calculate_reflection(ray) {
                Some(ans) => return Some(ans),
                None => 1,
            };
        }
        None
    }
}

/// Silly little enum for readability
pub enum DirChange {
    ChangeDir,
    KeepDir,
}

impl Ray {
    /// Initialise [`Ray`] struct from config vals
    pub fn init_from_cfgs(init_source: &RayInit, prog_config: &ProgConfig) -> Ray {
        let mut ray: Ray = Ray {
            range_vals: vec![0.0; prog_config.max_it + 1],
            depth_vals: vec![0.0; prog_config.max_it + 1],
            time_vals: vec![0.0; prog_config.max_it + 1],
            ray_param: init_source.init_ang.cos() / init_source.init_sound_speed,
            ray_iter: 0,
            ray_id: Uuid::new_v4().to_string(),
        };

        // initialise position and time
        ray.range_vals[0] = init_source.range_pos;
        ray.depth_vals[0] = init_source.depth_pos;
        ray.time_vals[0] = init_source.init_time;
        ray
    }

    /// Update step for ray
    pub fn update_iteration(
        &mut self,
        c_i: &f64,
        c_i1: &f64,
        g_i: &f64,
        depth_dir: &f64,
        depth_step: &f64,
    ) -> DirChange {
        if (self.ray_param * *c_i1).abs() < 1.0 {
            // calculate range and time steps for given depth step
            self.range_vals[self.ray_iter + 1] = self.range_vals[self.ray_iter]
                + (1.0 / (self.ray_param * g_i))
                    * ((1.0 - (self.ray_param * *c_i).powi(2)).sqrt()
                        - (1.0 - (self.ray_param * *c_i1).powi(2)).sqrt());
            self.time_vals[self.ray_iter + 1] = self.time_vals[self.ray_iter]
                + (((*c_i1 / *c_i) * (1.0 + (1.0 - (self.ray_param * *c_i).powi(2)).sqrt())
                    / (1.0 + (1.0 - (self.ray_param * c_i1).powi(2)).sqrt()))
                .ln()
                    / g_i)
                    .abs();
            // step depth values
            self.depth_vals[self.ray_iter + 1] =
                self.depth_vals[self.ray_iter] + *depth_dir * depth_step;
            DirChange::KeepDir
        } else {
            // Set depth to same value as ray is turning in this layer of medium
            self.depth_vals[self.ray_iter + 1] = self.depth_vals[self.ray_iter];
            // Calculate range and time steps for turning ray
            self.range_vals[self.ray_iter + 1] = self.range_vals[self.ray_iter]
                + 2.0 * (1.0 - (self.ray_param * *c_i).powi(2)).sqrt() / (self.ray_param * g_i);
            self.time_vals[self.ray_iter + 1] = self.time_vals[self.ray_iter]
                + 2.0
                    * ((1.0 + (1.0 - (self.ray_param * *c_i).powi(2)).sqrt())
                        / (self.ray_param * *c_i))
                        .ln()
                    / g_i.abs();
            DirChange::ChangeDir
        }
    }

    /// Remove unneeded cells
    pub fn truncate_ray(&mut self) {
        self.depth_vals.truncate(self.ray_iter + 1);
        self.range_vals.truncate(self.ray_iter + 1);
        self.time_vals.truncate(self.ray_iter + 1);
    }

    /// Update reflection step of ray propagation
    pub fn update_intersection(&mut self, reflect_ans: &ReflectResult, ssp: &Ssp) {
        // set next ray location to intersection location
        self.range_vals[self.ray_iter + 1] = reflect_ans.range;
        self.depth_vals[self.ray_iter + 1] = reflect_ans.depth;
        // reset xi ray parameter on new reflected source angle
        self.ray_param = reflect_ans.ang.cos() / ssp.interp_sound_speed(reflect_ans.depth);
        self.ray_iter += 1;
        // set next range and depth vals to minor offset from body to avoid ray getting
        // caught inside of [`Body`] struct
        self.range_vals[self.ray_iter + 1] =
            reflect_ans.range + REFLECT_OFFSET * reflect_ans.ang.cos();
        self.depth_vals[self.ray_iter + 1] =
            reflect_ans.depth + REFLECT_OFFSET * reflect_ans.ang.sin();
        self.time_vals[self.ray_iter + 1] = reflect_ans.time;
    }

    /// Trace ray using geometric theory
    pub fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        env_config: &EnvConfig,
        ssp: &Ssp,
    ) -> Self {
        let mut c_i: f64;
        let mut c_i1: f64;
        let mut g_i: f64;
        let ang: f64 = init_source.init_ang;
        let mut depth_dir: f64 = ang.sin().signum();
        // let mut range_dir: f64 = ang.cos().signum();
        let mut ray: Ray = Ray::init_from_cfgs(init_source, prog_config);
        // set initial sound speed values
        c_i = ssp.interp_sound_speed(ray.depth_vals[0]);

        for i in 0..env_config.isospaces.len() {
            if env_config.isospaces[i]
                .body
                .contains_point(&ray.range_vals[0], &ray.depth_vals[0])
            {
                ray.iso_trace(&env_config.isospaces[i]);
                ray.ray_iter += 1;
                ray.truncate_ray();
                return ray;
            }
        }

        while (ray.ray_iter < prog_config.max_it - init_source.init_iter)
            && (init_source
                .range_lims
                .contains(&ray.range_vals[ray.ray_iter]))
        {
            // set next SSP value
            c_i1 = ssp.interp_sound_speed(
                ray.depth_vals[ray.ray_iter] + depth_dir * prog_config.depth_step,
            );
            // calculate local sound speed gradient
            g_i = (c_i1 - c_i) / prog_config.depth_step;
            // iterate depth step. This function will update value of c_i and depth_dir as required
            match ray.update_iteration(&c_i, &c_i1, &g_i, &depth_dir, &prog_config.depth_step) {
                DirChange::KeepDir => c_i = c_i1,
                DirChange::ChangeDir => depth_dir = -depth_dir,
            };

            // Check for intersections on all bodies in simulation
            if let Some(reflect_ans) = env_config.check_all_body_reflections(&ray) {
                // Update intersection in step
                ray.update_intersection(&reflect_ans, ssp);
                // recalculate depth direction
                depth_dir = reflect_ans.ang.sin().signum();
                // Interpolate for next sound speed profile value
                c_i = ssp.interp_sound_speed(reflect_ans.depth);
            }
            // step iter value for calculation step taken
            ray.ray_iter += 1;
        }

        // truncate vectors to remove any wasted space
        ray.truncate_ray();
        ray
    }

    /// Trace straight ray through constant SSP areas of [`IsoSpace`]
    fn iso_trace(&mut self, iso_space: &IsoSpace) {
        //TODO: implement something to make max/min methods more terse
        let range_min: f64 = iso_space
            .body
            .range_vals
            .iter()
            .cloned()
            .fold(f64::INFINITY, |a, b| a.min(b));
        let range_max: f64 = iso_space
            .body
            .range_vals
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let depth_min: f64 = iso_space
            .body
            .depth_vals
            .iter()
            .cloned()
            .fold(f64::INFINITY, |a, b| a.min(b));
        let depth_max: f64 = iso_space
            .body
            .depth_vals
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);

        // Set ray endpoint outside of IsoSpace object
        self.range_vals[self.ray_iter + 1] = self.range_vals[self.ray_iter] + range_max - range_min;
        self.depth_vals[self.ray_iter + 1] = self.depth_vals[self.ray_iter] + depth_max - depth_min;

        // Update ray location and time with intersection
        match iso_space.body.check_finite_intersection(self) {
            Some(loc) => {
                self.range_vals[self.ray_iter + 1] = loc.range;
                self.depth_vals[self.ray_iter + 1] = loc.depth;
                self.time_vals[self.ray_iter + 1] =
                    ((self.range_vals[self.ray_iter + 1] - self.range_vals[self.ray_iter]).powi(2)
                        + (self.range_vals[self.ray_iter + 1] - self.range_vals[self.ray_iter])
                            .powi(2))
                    .sqrt()
                        / iso_space.sound_speed;
            }
            None => panic!(
                "No intersection found in IsoSpace object for ray {:}",
                self.ray_id
            ),
        }
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
        let ray_time_step: f64 = ray.time_vals[ray.ray_iter + 1] - ray.time_vals[ray.ray_iter];
        // Define valid ranges for solution parameters as defined in theory document
        let edge_param_range: RangeInclusive<f64> = RangeInclusive::new(0.0, 1.0);
        // iterate over each edge in polygon to find valid intersections
        for i in 0..(self.range_vals.len() - 1) {
            edge_depth_step = self.depth_vals[i + 1] - self.depth_vals[i];
            edge_range_step = self.range_vals[i + 1] - self.range_vals[i];
            // check for parallel edge and ray direction to save computation effort
            // TODO: Need else case
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
                    time: ray.time_vals[ray.ray_iter] + ray_dists_numeric[id] * ray_time_step,
                    edge_id: id,
                }),
                false => None,
            },
            None => None,
        }
    }

    fn calculate_reflection(&self, ray: &Ray) -> Option<ReflectResult> {
        let intersection_ans: IntersectLoc = self.check_finite_intersection(ray)?;
        let ray_ang: f64 = (ray.depth_vals[ray.ray_iter + 1] - ray.depth_vals[ray.ray_iter])
            .atan2(ray.range_vals[ray.ray_iter + 1] - ray.range_vals[ray.ray_iter]);
        let side_ang: f64 = self.get_edge_ang(intersection_ans.edge_id);
        let new_ang: f64 = 2.0 * side_ang - ray_ang;
        Some(ReflectResult {
            range: intersection_ans.range,
            depth: intersection_ans.depth,
            time: intersection_ans.time,
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

    /// check for intersection with side i
    fn intersects_dist(
        &self,
        ray_range_vals: [f64; 2],
        ray_depth_vals: [f64; 2],
        edge_id: usize,
    ) -> Option<f64> {
        let edge_param_range: RangeInclusive<f64> = RangeInclusive::new(0.0, 1.0);
        let edge_depth_vals: [f64; 2] = [self.depth_vals[edge_id], self.depth_vals[edge_id + 1]];
        let edge_range_vals: [f64; 2] = [self.range_vals[edge_id], self.range_vals[edge_id + 1]];
        let edge_depth_step: f64 = edge_depth_vals[1] - edge_depth_vals[0];
        let edge_range_step: f64 = edge_range_vals[1] - edge_range_vals[0];
        let ray_range_step: f64 = ray_range_vals[1] - ray_range_vals[0];
        let ray_depth_step: f64 = ray_depth_vals[1] - ray_depth_vals[0];
        if (edge_range_step / edge_depth_step) != (ray_range_step / ray_depth_step) {
            let edge_dist: f64 = ((self.range_vals[edge_id] - ray_range_vals[0]) * ray_depth_step
                - (self.depth_vals[edge_id] - ray_depth_vals[0]) * ray_range_step)
                / ((self.depth_vals[edge_id + 1] - self.depth_vals[edge_id]) * ray_range_step
                    - (self.range_vals[edge_id + 1] - self.range_vals[edge_id]) * ray_depth_step);
            if edge_param_range.contains(&edge_dist) {
                return match ray_range_step != 0.0 {
                    true => intersection_ray_dist_param(
                        edge_dist,
                        edge_range_step,
                        ray_range_step,
                        edge_range_vals[0],
                        ray_range_vals[0],
                    ),
                    false => intersection_ray_dist_param(
                        edge_dist,
                        edge_depth_step,
                        ray_depth_step,
                        edge_range_vals[0],
                        ray_depth_vals[0],
                    ),
                };
            }
        }
        None
    }

    /// Calculates if point at range [m] and depth [m] is inside of polygon shape
    fn contains_point(&self, range: &f64, depth: &f64) -> bool {
        // cast vertical ray above surface
        let ray_range_vals: [f64; 2] = [*range, *range];
        let ray_depth_vals: [f64; 2] = [*depth, -1.0];
        let mut valid_count: usize = 0;
        for i in 0..(self.range_vals.len() - 1) {
            if self
                .intersects_dist(ray_range_vals, ray_depth_vals, i)
                .is_some()
            {
                valid_count += 1;
            }
        }
        (valid_count % 2) == 1
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
pub fn trace_rays(cfg: Config) -> Vec<Ray> {
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
        .map(|ray_source| {
            Ray::trace_from_init_source(
                ray_source,
                &cfg.prog_config,
                &cfg.env_config,
                &cfg.env_config.ssp,
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
