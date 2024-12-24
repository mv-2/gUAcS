use crate::path_tracing::{Body, HalfSpace, Ssp};
use serde::Deserialize;
use std::fs::File;

pub const OUTPUT_DIR: &str = "output_data";
pub const CONFIG_PATH: &str = "configs";

/// Overall config opject to aid in loading serialized config jsons
#[derive(Deserialize, Debug)]
pub struct Config {
    pub prog_config: ProgConfig,
    pub env_config: EnvConfig,
    pub sources: Vec<SourceConfig>,
}

/// Config to store programmatic data not relevant to theory of simulation
#[derive(Deserialize, Debug)]
pub struct ProgConfig {
    pub max_it: usize,
    pub depth_step: f64,
    pub save_to_csv: bool,
    pub max_range: f64,
}

/// Stores environmental constant data for simulation (SSP and density profile information)
#[derive(Deserialize, Debug)]
pub struct EnvConfig {
    pub ssp: Ssp,
    pub swell_height: f64,
    pub bodies: Vec<Body>,
    pub halfspaces: Vec<HalfSpace>,
}

/// Stores information of single source in sound field
#[derive(Deserialize, Debug)]
pub struct SourceConfig {
    // eventually each ray should be set a source_level based on its angle and the type of shot
    pub range_pos: f64,
    pub depth_pos: f64,
    pub ray_fan_limits: [f64; 2],
    pub n_rays: usize,
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
