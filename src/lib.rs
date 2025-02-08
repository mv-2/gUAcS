use pyo3::prelude::*;
pub mod interface;
pub mod math_util;
pub mod path_tracing;

use crate::interface::Config;
use crate::path_tracing::{trace_from_config, Ray};
use glob::glob;
use std::fs::remove_file;

#[pyfunction]
#[pyo3(name = "run_sim")]
fn python_rays(config_path: String, output_path: String) -> PyResult<Vec<Ray>> {
    // TODO: Possibly remove this later on as deleting the users files is suboptimal
    for path in glob(&format!("{output_path}/*.csv")[..]).expect("Failed to find output_data files")
    {
        match path {
            Ok(path) => remove_file(path).unwrap(),
            Err(e) => println!("{:?}", e),
        }
    }
    let cfg: Config = Config::from_json(config_path);
    Ok(trace_from_config(cfg, output_path))
}

#[pymodule]
fn guacs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(python_rays, m)?)?;
    Ok(())
}
