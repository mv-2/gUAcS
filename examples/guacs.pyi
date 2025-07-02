from typing import List, Tuple

def trace_rays(config: Config) -> List[Ray]: ...
def trace_beams(config: Config) -> List[PyBeam]: ...

class Config:
    prog_config: ProgConfig
    env_config: EnvConfig
    sources: List[SourceConfig]

class ProgConfig:
    max_it: int
    depth_step: float
    max_range: float
    min_range: float
    output_path: str
    pq_solver: str

class EnvConfig:
    ssp: Ssp
    swell_height: float
    bodies: List[Body]
    isospaces: List[IsoSpace]

class Ssp:
    ssp_knots: List[float]
    ssp_coefs: List[float]
    ssp_degree: int

class Body:
    range_vals: List[float]
    depth_vals: List[float]

class IsoSpace:
    body: Body
    sound_speed: float
    density: float

class SourceConfig:
    range_pos: float
    depth_pos: float
    ray_fan_limits: Tuple[float, float]
    n_rays: int
    source_level: float
    frequency: float

class Ray:
    range_vals: List[float]
    depth_vals: List[float]
    time_vals: List[float]
    ray_param: float
    ray_iter: int
    ray_id: str

class PyBeam:
    central_ray: Ray
    p_re: List[float]
    q_re: List[float]
    p_im: List[float]
    q_im: List[float]
