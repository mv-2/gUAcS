use crate::beam_tracing::{PressureField, PyBeam};
use crate::path_tracing::{Body, Ssp};

use pyo3::prelude::*;
use serde::Deserialize;

/// Overall config opject to aid in loading serialized config jsons
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct RayConfig {
    #[pyo3(get, set)]
    pub prog_config: ProgConfig,
    #[pyo3(get, set)]
    pub env_config: EnvConfig,
    #[pyo3(get, set)]
    pub sources: Vec<SourceConfig>,
}

/// Beam Config struct for python
#[derive(Clone)]
#[pyclass]
pub struct BeamConfig {
    #[pyo3(get, set)]
    ray_config: RayConfig,
    #[pyo3(get, set)]
    pq_solver: String,
    #[pyo3(get, set)]
    pressure_locs: Vec<(f64, f64)>,
    #[pyo3(get, set)]
    window_width: Option<f64>,
}

/// Beam Config struct for rust
#[derive(Clone)]
pub struct BeamConfigRust {
    pub ray_config: RayConfig,
    pub pq_solver: SolverMethod,
    pub pressure_locs: Vec<(f64, f64)>,
    pub window_width: Option<f64>,
}

impl From<BeamConfig> for BeamConfigRust {
    fn from(beam_config: BeamConfig) -> Self {
        BeamConfigRust {
            ray_config: beam_config.ray_config,
            pq_solver: SolverMethod::from(beam_config.pq_solver),
            pressure_locs: beam_config.pressure_locs,
            window_width: beam_config.window_width,
        }
    }
}

/// RayConfig to store programmatic data not relevant to theory of simulation
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
    #[pyo3(get, set)]
    pub isospaces: Vec<IsoSpace>,
}

/// Stores information of single source in sound field
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct SourceConfig {
    // eventually each ray should be set a source_level based on its angle and the type of shot
    #[pyo3(get, set)]
    pub range_pos: f64,
    #[pyo3(get, set)]
    pub depth_pos: f64,
    #[pyo3(get, set)]
    pub ray_fan_limits: [f64; 2],
    #[pyo3(get, set)]
    pub n_rays: usize,
    #[pyo3(get, set)]
    pub source_level: f64,
    #[pyo3(get, set)]
    pub frequency: f64,
}

/// Constant property space. Implementation will be usable for RF simulation as well
#[derive(Deserialize, Debug, Clone)]
#[pyclass]
pub struct IsoSpace {
    #[pyo3(get, set)]
    pub body: Body,
    #[pyo3(get, set)]
    pub sound_speed: f64,
    #[pyo3(get, set)]
    pub density: f64,
}

/// Python equivalent to [`PressureField`] object to allow for Complex values to be transmitted
/// correctly
#[derive(Clone)]
#[pyclass]
pub struct PressureFieldPy {
    #[pyo3(get, set)]
    pub locations: Vec<(f64, f64)>,
    #[pyo3(get, set)]
    pub re: Vec<f64>,
    #[pyo3(get, set)]
    pub im: Vec<f64>,
}

/// Make py object from rust struct
impl From<PressureField> for PressureFieldPy {
    fn from(pressure_field: PressureField) -> Self {
        PressureFieldPy {
            locations: pressure_field.locations,
            re: pressure_field.pressures.iter().map(|&z| z.re).collect(),
            im: pressure_field.pressures.iter().map(|&z| z.im).collect(),
        }
    }
}

#[pyclass]
pub struct BeamResult {
    #[pyo3(get, set)]
    pub beams: Vec<PyBeam>,
    #[pyo3(get, set)]
    pub pressures: PressureFieldPy,
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
        frequency: f64,
    ) -> Self {
        SourceConfig {
            range_pos,
            depth_pos,
            ray_fan_limits,
            n_rays,
            source_level,
            frequency,
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
impl IsoSpace {
    #[new]
    fn py_new(body: Body, sound_speed: f64, density: f64) -> Self {
        IsoSpace {
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
    fn py_new(bodies: Vec<Body>, ssp: Ssp, swell_height: f64, isospaces: Vec<IsoSpace>) -> Self {
        EnvConfig {
            bodies,
            ssp,
            swell_height,
            isospaces,
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
impl RayConfig {
    #[new]
    fn py_new(env_config: EnvConfig, prog_config: ProgConfig, sources: Vec<SourceConfig>) -> Self {
        RayConfig {
            env_config,
            prog_config,
            sources,
        }
    }
}

#[pymethods]
impl BeamConfig {
    #[new]
    fn py_new(
        ray_config: RayConfig,
        pq_solver: String,
        pressure_locs: Vec<(f64, f64)>,
        window_width: Option<f64>,
    ) -> Self {
        BeamConfig {
            ray_config,
            pq_solver,
            pressure_locs,
            window_width,
        }
    }
}

#[pymethods]
impl PressureFieldPy {
    #[new]
    fn py_new(locations: Vec<(f64, f64)>, re: Vec<f64>, im: Vec<f64>) -> Self {
        PressureFieldPy { locations, re, im }
    }
}

#[pymethods]
impl BeamResult {
    #[new]
    fn py_new(beams: Vec<PyBeam>, pressures: PressureFieldPy) -> Self {
        BeamResult { beams, pressures }
    }
}

// Enum to set solver method for p-q equations
#[derive(Clone, Copy)]
pub enum SolverMethod {
    RungeKutta4,
    Radau3IA,
    BackEuler,
    Radau3IIA,
}

impl From<String> for SolverMethod {
    fn from(s: String) -> Self {
        match s.as_str() {
            "BackwardEuler" => SolverMethod::BackEuler,
            "RungeKutta4" => SolverMethod::RungeKutta4,
            "Radau3IA" => SolverMethod::Radau3IA,
            "Radau3IIA" => SolverMethod::Radau3IIA,
            _ => panic!("Invalid p-q solver method: {s}"),
        }
    }
}
