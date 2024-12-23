from pydantic.dataclasses import dataclass
from pydantic import BaseModel
from typing import List, Tuple
import json
import numpy as np


class Ray:
    def __init__(self, path: str):
        ray_id = path.split("/")[-1].split(".csv")[0]
        data = np.loadtxt(path, delimiter=",", dtype=float)
        self.range_vals = data[:, 0]
        self.depth_vals = data[:, 1]
        self.time_vals = data[:, 2]
        self.ray_id = ray_id


@dataclass
class Body:
    range_vals: List[float]
    depth_vals: List[float]


@dataclass
class HalfSpace:
    body: Body
    sound_speed: float
    density: float


@dataclass
class SourceConfig:
    range_pos: float
    depth_pos: float
    ray_fan_limits: Tuple[float, float]
    n_rays: int
    source_level: float


@dataclass
class ProgConfig:
    max_it: int
    depth_step: float
    save_to_csv: bool


@dataclass
class Ssp:
    ssp_knots: list
    ssp_coefs: list
    ssp_degree: int


@dataclass
class EnvConfig:
    ssp: Ssp
    swell_height: float
    bodies: List[Body]
    halfspaces: List[HalfSpace]


class Config(BaseModel):
    prog_config: ProgConfig
    env_config: EnvConfig
    sources: List[SourceConfig]

    def write_to_json(self, save_path: str):
        with open(save_path, "w") as f:
            json.dump(json.loads(self.json()), f)
