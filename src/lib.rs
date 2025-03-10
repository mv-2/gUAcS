use path_tracing::HalfSpace;
use pyo3::prelude::*;
pub mod interface;
pub mod math_util;
pub mod path_tracing;

use crate::interface::*;
use crate::path_tracing::{trace_from_config, Body, Ray, Ssp};
use glob::glob;
use std::fs::remove_file;

#[pyfunction]
#[pyo3(name = "run_sim")]
fn python_rays(config: Config, output_path: String) -> PyResult<Vec<Ray>> {
    // TODO: Possibly remove this later on as deleting the users files is suboptimal
    for path in glob(&format!("{output_path}/*.csv")[..]).expect("Failed to find output_data files")
    {
        match path {
            Ok(path) => remove_file(path).unwrap(),
            Err(e) => println!("{:?}", e),
        }
    }
    Ok(trace_from_config(config, output_path))
}

#[pymodule]
fn guacs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(python_rays, m)?)?;
    m.add_class::<Config>()?;
    m.add_class::<ProgConfig>()?;
    m.add_class::<EnvConfig>()?;
    m.add_class::<SourceConfig>()?;
    m.add_class::<Ssp>()?;
    m.add_class::<Body>()?;
    m.add_class::<HalfSpace>()?;
    m.add_class::<Ray>()?;
    Ok(())
}
