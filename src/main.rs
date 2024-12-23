pub mod interface;
pub mod math_util;
pub mod path_tracing;

use crate::interface::{Config, CONFIG_PATH};
use crate::path_tracing::{trace_from_config, Ray};
use glob::glob;
use std::fs::remove_file;

fn main() {
    for path in glob("output_data/rays/*.csv").expect("Failed to find output_data files") {
        match path {
            Ok(path) => remove_file(path).unwrap(),
            Err(e) => println!("{:?}", e),
        }
    }
    let cfg: Config = Config::from_json(format!("{CONFIG_PATH}/run_config.json"));
    let _: Vec<Ray> = trace_from_config(cfg);
}
