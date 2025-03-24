use crate::path_tracing::{Body, HalfSpace, Ssp};
use pyo3::prelude::*;
use serde::Deserialize;
use std::fs::File;

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
    pub save_to_csv: bool,
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
    #[pyo3(get, set)]
    pub halfspaces: Vec<HalfSpace>,
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

impl Config {
    /// unpack json config into nested rust structs
    pub fn from_json(file_path: String) -> Self {
        // This kind of sucks due to all of the cloning but you can't argue with results
        let file: File = File::open(file_path).expect("File should open read only");
        let json: serde_json::Value =
            serde_json::from_reader(file).expect("File should be well formed json");

        let env_val: serde_json::Value = json
            .get("env_config")
            .expect("Error when reading env_config field")
            .clone();

        let prog_val: serde_json::Value = json
            .get("prog_config")
            .expect("Error when reading prog_config field")
            .clone();

        let sources: serde_json::Value = json
            .get("sources")
            .expect("Error when reading sources field")
            .clone();

        let prog_config: ProgConfig = match serde_json::from_value(prog_val) {
            Ok(cfg) => cfg,
            Err(_) => panic!("Failed deserializing prog_config"),
        };

        let ssp_val: serde_json::Value = env_val
            .get("ssp")
            .expect("Error when reading SSP field")
            .clone();

        let ssp: Ssp = match serde_json::from_value(ssp_val) {
            Ok(ssp) => ssp,
            Err(_) => panic!("Failed deserializing ssp"),
        };

        let swell_val: serde_json::Value = env_val
            .get("swell_height")
            .expect("Error when reading swell_height")
            .clone();

        let swell_height: f64 = match serde_json::from_value(swell_val) {
            Ok(swell) => swell,
            Err(_) => panic!("Failed deserializing swell_height"),
        };

        let body_arr_val: serde_json::Value = env_val
            .get("bodies")
            .expect("Error when reading body array")
            .clone();

        let body_vec: Vec<Body> = match serde_json::from_value(body_arr_val) {
            Ok(bds) => bds,
            Err(_) => panic!("Error when deserializing body array"),
        };

        let halfspace_arr_val: serde_json::Value = env_val
            .get("halfspaces")
            .expect("Error when reading body array")
            .clone();

        let halfspace_vec: Vec<HalfSpace> = match serde_json::from_value(halfspace_arr_val) {
            Ok(hfsp) => hfsp,
            Err(_) => panic!("Error when deserializing halfspace array"),
        };

        let env_config: EnvConfig = EnvConfig {
            ssp,
            swell_height,
            bodies: body_vec,
            halfspaces: halfspace_vec,
        };

        let source_vec: Vec<SourceConfig> = match serde_json::from_value(sources) {
            Ok(sces) => sces,
            Err(_) => panic!("Error when deserializing source array"),
        };

        Config {
            prog_config,
            env_config,
            sources: source_vec,
        }
    }
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

#[pymethods]
impl HalfSpace {
    #[new]
    fn py_new(body: Body, sound_speed: f64, density: f64) -> Self {
        HalfSpace {
            body,
            sound_speed,
            density,
        }
    }
}

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
    fn py_new(bodies: Vec<Body>, ssp: Ssp, swell_height: f64, halfspaces: Vec<HalfSpace>) -> Self {
        EnvConfig {
            bodies,
            ssp,
            swell_height,
            halfspaces,
        }
    }
}

#[pymethods]
impl ProgConfig {
    #[new]
    fn py_new(
        depth_step: f64,
        max_it: usize,
        save_to_csv: bool,
        max_range: f64,
        min_range: f64,
        output_path: String,
    ) -> Self {
        ProgConfig {
            depth_step,
            max_it,
            save_to_csv,
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
