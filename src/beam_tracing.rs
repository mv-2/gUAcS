use crate::interface::{EnvConfig, ProgConfig};
use crate::path_tracing::{Ray, RayInit, Ssp};
use num::complex::Complex64;

// Stores Beam propagation
// #[pyclass]
// #[derive(Debug)]
pub struct Beam {
    // #[pyo3(get, set)]
    pub central_ray: Ray,
    // #[pyo3(get, set)]
    pub p: Complex64,
    // #[pyo3(get, set)]
    pub q: Complex64,
}

impl Beam {
    fn init_from_configs(init_source: &RayInit, prog_config: &ProgConfig) -> Beam {
        Beam {
            central_ray: Ray::init_from_cfgs(init_source, prog_config),
            p: Complex64::new(1.0, 0.0),
            q: Complex64::new(0.0, 1.0 / init_source.init_sound_speed),
        }
    }

    fn trace_from_init_source(
        init_source: &RayInit,
        prog_config: &ProgConfig,
        env_config: &EnvConfig,
        ssp: &Ssp,
    ) {
        let beam: Beam = Beam::init_from_configs(init_source, prog_config);
        todo!();
    }
}
