use crate::path_tracing::{Body, Ssp};
use pyo3::prelude::*;
use serde::Deserialize;

/// Overall config opject to aid in loading serialized config jsons
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct Config {
    #[pyo3(get, set)]
    pub prog_config: ProgConfig,
    #[pyo3(get, set)]
    pub env_config: EnvConfig,
    #[pyo3(get, set)]
    pub sources: Vec<SourceConfig>,
}

/// Config to store programmatic data not relevant to theory of simulation
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct ProgConfig {
    #[pyo3(get, set)]
    pub max_it: usize,
    #[pyo3(get, set)]
    pub depth_step: f64,
    #[pyo3(get, set)]
    pub max_range: f64,
    #[pyo3(get, set)]
    pub min_range: f64,
    #[pyo3(get, set)]
    pub output_path: String,
}

/// Stores environmental constant data for simulation (SSP and density profile information)
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct EnvConfig {
    #[pyo3(get, set)]
    pub ssp: Ssp,
    #[pyo3(get, set)]
    pub swell_height: f64,
    #[pyo3(get, set)]
    pub bodies: Vec<Body>,
    // #[pyo3(get, set)]
    // pub halfspaces: Vec<HalfSpace>,
}

/// Stores information of single source in sound field
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct SourceConfig {
    // eventually each ray should be set a source_level based on its angle and the type of shot
    // #[pyo3(get, set)]
    pub range_pos: f64,
    #[pyo3(get, set)]
    pub depth_pos: f64,
    #[pyo3(get, set)]
    pub ray_fan_limits: [f64; 2],
    #[pyo3(get, set)]
    pub n_rays: usize,
    #[pyo3(get, set)]
    pub source_level: f64,
}

// Python __new__ constructors for required structs
#[pymethods]
impl SourceConfig {
    #[new]
    fn py_new(
        range_pos: f64,
        depth_pos: f64,
        ray_fan_limits: [f64; 2],
        n_rays: usize,
        source_level: f64,
    ) -> Self {
        SourceConfig {
            range_pos,
            depth_pos,
            ray_fan_limits,
            n_rays,
            source_level,
        }
    }
}

#[pymethods]
impl Body {
    #[new]
    fn py_new(range_vals: Vec<f64>, depth_vals: Vec<f64>) -> Self {
        Body {
            range_vals,
            depth_vals,
        }
    }
}

// #[pymethods]
// impl HalfSpace {
//     #[new]
//     fn py_new(body: Body, sound_speed: f64, density: f64) -> Self {
//         HalfSpace {
//             body,
//             sound_speed,
//             density,
//         }
//     }
// }

#[pymethods]
impl Ssp {
    #[new]
    fn py_new(ssp_knots: Vec<f64>, ssp_coefs: Vec<f64>, ssp_degree: usize) -> Self {
        Ssp {
            ssp_knots,
            ssp_coefs,
            ssp_degree,
        }
    }
}

#[pymethods]
impl EnvConfig {
    #[new]
    fn py_new(
        bodies: Vec<Body>,
        ssp: Ssp,
        swell_height: f64,
        //halfspaces: Vec<HalfSpace>
    ) -> Self {
        EnvConfig {
            bodies,
            ssp,
            swell_height,
            // halfspaces,
        }
    }
}

#[pymethods]
impl ProgConfig {
    #[new]
    fn py_new(
        depth_step: f64,
        max_it: usize,
        max_range: f64,
        min_range: f64,
        output_path: String,
    ) -> Self {
        ProgConfig {
            depth_step,
            max_it,
            max_range,
            min_range,
            output_path,
        }
    }
}

#[pymethods]
impl Config {
    #[new]
    fn py_new(env_config: EnvConfig, prog_config: ProgConfig, sources: Vec<SourceConfig>) -> Self {
        Config {
            env_config,
            prog_config,
            sources,
        }
    }
}
