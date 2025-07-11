use crate::interface::{EnvConfig, IsoSpace, ProgConfig, RayConfig, SourceConfig};
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
    #[pyo3(get, set)]
    pub frequency: f64,
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
    pub ang_step: f64,
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

        let ang_step: f64 = match source.n_rays == 1 {
            true => 2.0 * f64::consts::PI,
            false => (source.ray_fan_limits[1] - source.ray_fan_limits[0]) / source.n_rays as f64,
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
            ang_step,
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

/// Stores result of reflection calculation in [`Body::calculate_reflection()`]
#[derive(Debug)]
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
        // I would love to directly iterate body values but unfortunately borrow checker is whiny
        let test_range_vals: [f64; 2] = [
            ray.range_vals[ray.ray_iter],
            ray.range_vals[ray.ray_iter + 1],
        ];
        let test_depth_vals: [f64; 2] = [
            ray.depth_vals[ray.ray_iter],
            ray.depth_vals[ray.ray_iter + 1],
        ];
        let test_time_vals: [f64; 2] =
            [ray.time_vals[ray.ray_iter], ray.time_vals[ray.ray_iter + 1]];
        let ray_ang: f64 =
            (ray.ray_param * self.ssp.interp_sound_speed(ray.depth_vals[ray.ray_iter])).acos();
        for i in 0..self.bodies.len() {
            match self.bodies[i].calculate_reflection(
                &test_range_vals,
                &test_depth_vals,
                &test_time_vals,
                &ray_ang,
            ) {
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
            frequency: init_source.frequency,
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
                ray.iso_trace(&env_config.isospaces[i], prog_config.max_it, &ang);
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
    pub fn iso_trace(&mut self, iso_space: &IsoSpace, max_it: usize, ang: &f64) {
        // this is dumb
        let mut ang: f64 = *ang;
        let mut test_range_vals: [f64; 2];
        let mut test_depth_vals: [f64; 2];
        let mut test_time_vals: [f64; 2];
        //TODO: implement something to make max/min methods more terse
        let range_max: f64 = iso_space
            .body
            .range_vals
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| f64::max(a, b));
        let range_min: f64 = iso_space
            .body
            .range_vals
            .iter()
            .fold(f64::INFINITY, |a, &b| f64::min(a, b));
        let depth_max: f64 = iso_space
            .body
            .depth_vals
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| f64::max(a, b));
        let depth_min: f64 = iso_space
            .body
            .depth_vals
            .iter()
            .fold(f64::INFINITY, |a, &b| f64::min(a, b));
        let test_dist: f64 =
            1.1 * ((range_max - range_min).powi(2) + (depth_max - depth_min).powi(2)).sqrt();
        while self.ray_iter < max_it {
            // Set ray endpoint outside of IsoSpace obj
            test_range_vals = [
                self.range_vals[self.ray_iter] + REFLECT_OFFSET * ang.cos(),
                self.range_vals[self.ray_iter] + test_dist * ang.cos(),
            ];
            test_depth_vals = [
                self.depth_vals[self.ray_iter] + REFLECT_OFFSET * ang.sin(),
                self.depth_vals[self.ray_iter] + test_dist * ang.sin(),
            ];
            test_time_vals = [self.time_vals[self.ray_iter], 0_f64];

            // Update ray location and time with intersection
            match iso_space.body.calculate_reflection(
                &test_range_vals,
                &test_depth_vals,
                &test_time_vals,
                &ang,
            ) {
                Some(loc) => {
                    self.range_vals[self.ray_iter + 1] = loc.range;
                    self.depth_vals[self.ray_iter + 1] = loc.depth;
                    self.time_vals[self.ray_iter + 1] = self.time_vals[self.ray_iter]
                        + ((self.range_vals[self.ray_iter + 1] - self.range_vals[self.ray_iter])
                            .powi(2)
                            + (self.range_vals[self.ray_iter + 1]
                                - self.range_vals[self.ray_iter])
                                .powi(2))
                        .sqrt()
                            / iso_space.sound_speed;
                    ang = loc.ang;
                }
                None => panic!("No intesection found in IsoSpace for ray {:}", self.ray_id),
            }
            self.ray_iter += 1;
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
    fn check_finite_intersection(
        &self,
        range_vals: &[f64],
        depth_vals: &[f64],
        time_vals: &[f64],
    ) -> Option<IntersectLoc> {
        let mut edge_dist: f64;
        let mut edge_depth_step: f64;
        let mut edge_range_step: f64;
        let mut ray_dist_vals: Vec<Option<f64>> = vec![None; self.range_vals.len() - 1];
        let ray_range_step: f64 = range_vals[1] - range_vals[0];
        let ray_depth_step: f64 = depth_vals[1] - depth_vals[0];
        let ray_time_step: f64 = time_vals[1] - time_vals[0];
        // Define valid ranges for solution parameters as defined in theory document
        let edge_param_range: RangeInclusive<f64> = RangeInclusive::new(0.0, 1.0);
        // iterate over each edge in polygon to find valid intersections
        for i in 0..(self.range_vals.len() - 1) {
            edge_depth_step = self.depth_vals[i + 1] - self.depth_vals[i];
            edge_range_step = self.range_vals[i + 1] - self.range_vals[i];
            // check for parallel edge and ray direction to save computation effort
            // TODO: Need else case
            if edge_range_step / edge_depth_step != ray_range_step / ray_depth_step {
                edge_dist = ((self.range_vals[i] - range_vals[0]) * ray_depth_step
                    - (self.depth_vals[i] - depth_vals[0]) * ray_range_step)
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
                            range_vals[0],
                        ),
                        false => intersection_ray_dist_param(
                            edge_dist,
                            edge_depth_step,
                            ray_depth_step,
                            self.depth_vals[i],
                            depth_vals[0],
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
                    range: range_vals[0] + ray_dists_numeric[id] * ray_range_step,
                    depth: depth_vals[0] + ray_dists_numeric[id] * ray_depth_step,
                    time: time_vals[0] + ray_dists_numeric[id] * ray_time_step,
                    edge_id: id,
                }),
                false => None,
            },
            None => None,
        }
    }

    pub fn calculate_reflection(
        &self,
        range_vals: &[f64],
        depth_vals: &[f64],
        time_vals: &[f64],
        ray_ang: &f64,
    ) -> Option<ReflectResult> {
        let intersection_ans: IntersectLoc =
            self.check_finite_intersection(range_vals, depth_vals, time_vals)?;
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

    /// Calculates if point at range and depth (in metres) is inside of polygon shape
    pub fn contains_point(&self, range: &f64, depth: &f64) -> bool {
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

/// Driver function to trace all rays from [`RayConfig`] struct
pub fn trace_rays(cfg: RayConfig) -> Vec<Ray> {
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
            frequency: 0.0,
        };
        let ans: IntersectLoc = test_body
            .check_finite_intersection(
                &test_ray.range_vals,
                &test_ray.depth_vals,
                &test_ray.time_vals,
            )
            .expect("No solution found");
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
            frequency: 0.0,
        };
        let ans: Option<IntersectLoc> = test_body.check_finite_intersection(
            &test_ray.range_vals,
            &test_ray.depth_vals,
            &test_ray.time_vals,
        );
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
            frequency: 0.0,
        };
        let ans: Option<IntersectLoc> = test_body.check_finite_intersection(
            &test_ray.range_vals,
            &test_ray.depth_vals,
            &test_ray.time_vals,
        );
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
            frequency: 0.0,
        };
        let ans: Option<IntersectLoc> = test_body.check_finite_intersection(
            &test_ray.range_vals,
            &test_ray.depth_vals,
            &test_ray.time_vals,
        );
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
            frequency: 0.0,
        };
        // TODO: VERIFY THIS
        let ray_ang: f64 = (test_ray.range_vals[1] - test_ray.range_vals[0])
            .atan2(test_ray.depth_vals[1] - test_ray.depth_vals[0]);
        let ang_out: ReflectResult = test_body
            .calculate_reflection(
                &test_ray.range_vals,
                &test_ray.depth_vals,
                &test_ray.time_vals,
                &ray_ang,
            )
            .unwrap();
        println!("{:}", ang_out.ang);
    }
}
