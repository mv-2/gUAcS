pub mod interface;
pub mod math_util;
pub mod path_tracing;

use crate::interface::Config;
use crate::path_tracing::{trace_from_config, Ray};

const CONFIG_PATH: &str = "configs";

fn main() {
    let cfg: Config = Config::from_json(format!("{CONFIG_PATH}/run_config.json"));
    let rays: Vec<Ray> = trace_from_config(cfg);
    println!("{:#?}", rays[0]);
}
