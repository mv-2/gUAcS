pub mod beam_tracing;
pub mod interface;
pub mod math_util;
pub mod path_tracing;

use beam_tracing::evaluate_field;
use pyo3::prelude::*;

use crate::beam_tracing::{trace_beams, Beam, BeamDisplay, PressureField, PyBeam};
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
fn python_beams(config: BeamConfig) -> PyResult<BeamResult> {
    let rust_config: BeamConfigRust = BeamConfigRust::from(config);
    // config 2 electric boogaloo
    let body_clone: Vec<Body> = rust_config.ray_config.env_config.bodies.clone();
    let pressure_locs: Vec<(f64, f64)> = rust_config.pressure_locs.clone();
    let ssp: Ssp = rust_config.ray_config.env_config.ssp.clone();
    let beams: Vec<Beam> = trace_beams(rust_config.clone());
    let pressures: PressureField = evaluate_field(rust_config, &beams);

    // Egui
    let beam_display: BeamDisplay = BeamDisplay {
        pressures: (0..pressures.pressures.len())
            .map(|i| pressures.pressures[i].norm())
            .collect(),
        ssp,
        ssp_pts: (0..5000).map(f64::from).collect(),
        cols: pressure_locs
            .iter()
            .filter(|x| x.0 == pressure_locs[0].0)
            .count(),
        rays: beams.iter().map(|bm| bm.central_ray.clone()).collect(),
        bodies: body_clone,
    };

    let native_options = eframe::NativeOptions::default();

    let _ = eframe::run_native(
        BeamDisplay::name(),
        native_options,
        Box::new(|_| Ok(Box::new(beam_display))),
    );

    Ok(BeamResult {
        // TEST: check if par_iter() is worth the call here (probably not)
        beams: beams.iter().map(|bm| PyBeam::from_beam(bm)).collect(),
        pressures: PressureFieldPy::from(pressures),
    })
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
    m.add_class::<BeamResult>()?;
    m.add_class::<PressureFieldPy>()?;
    Ok(())
}
