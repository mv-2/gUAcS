pub mod beam_tracing;
pub mod interface;
pub mod math_util;
pub mod path_tracing;

use pyo3::prelude::*;
use rayon::prelude::*;

use crate::beam_tracing::{trace_beams, PyBeam};
use crate::interface::*;
use crate::path_tracing::{trace_rays, Body, Ray, Ssp};

#[pyfunction]
#[pyo3(name = "trace_rays")]
fn python_rays(config: RayConfig) -> PyResult<Vec<Ray>> {
    // silly wrapper for python stuff
    Ok(trace_rays(config))
}

#[pyfunction]
#[pyo3(name = "trace_beams")]
#[allow(clippy::redundant_closure)]
fn python_beams(config: BeamConfig) -> PyResult<Vec<PyBeam>> {
    let rust_config: BeamConfigRust = BeamConfigRust::from(config);
    Ok(trace_beams(rust_config)
        .par_iter()
        .map(|bm| PyBeam::from_beam(bm))
        .collect())
}

#[pymodule]
fn guacs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(python_rays, m)?)?;
    m.add_function(wrap_pyfunction!(python_beams, m)?)?;
    m.add_class::<RayConfig>()?;
    m.add_class::<BeamConfig>()?;
    m.add_class::<ProgConfig>()?;
    m.add_class::<EnvConfig>()?;
    m.add_class::<SourceConfig>()?;
    m.add_class::<Ssp>()?;
    m.add_class::<Body>()?;
    m.add_class::<IsoSpace>()?;
    m.add_class::<Ray>()?;
    m.add_class::<PyBeam>()?;
    Ok(())
}
