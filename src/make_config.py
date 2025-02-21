from python_interface_tools import (
    Body,
    HalfSpace,
    SourceConfig,
    EnvConfig,
    Config,
    Ssp,
    ProgConfig,
)
from scipy.interpolate import splrep
import numpy as np


def munk_profile(depth: float) -> float:
    return 1500.0 * (
        1.0
        + 0.00737 * (2 * (depth - 1300) / 1500 - 1 + np.exp(-2 * (depth - 1300) / 1500))
    )


CONFIG_PATH = "configs"


if __name__ == "__main__":
    sources = [
        SourceConfig(
            range_pos=0.0,
            depth_pos=1000.0,
            ray_fan_limits=(-0.1, 0.1),
            n_rays=2,
            source_level=150,
        )
    ]

    bodies = [
        Body(
            range_vals=[20000.0, 20100.0, 20100.0, 20000.0, 20000.0],
            depth_vals=[1300.0, 1300.0, 1310.0, 1310.0, 1300.0],
        ),
        Body(
            range_vals=[0.0, 200000.0, 200000.0, 0.0, 0.0],
            depth_vals=[5000.0, 5000.0, 5100.0, 5100.0, 5000.0],
        ),
        Body(
            range_vals=[0.0, 200000.0, 200000.0, 0.0, 0.0],
            depth_vals=[0.0, 0.0, -1.0, -1.0, 0.0],
        ),
    ]

    halfspaces = [
        HalfSpace(
            body=Body(
                range_vals=[0.0, 100000.0, 100000.0, 0.0, 0.0],
                depth_vals=[5000.0, 5000.0, 5100.0, 5100.0, 5000.0],
            ),
            sound_speed=1600.0,
            density=1.5,
        )
    ]

    munk_depths = range(0, 5000, 10)
    munk_vals = [munk_profile(z) for z in munk_depths]

    (ssp_knots, ssp_coefs, ssp_degree) = splrep(munk_depths, munk_vals, k=3)

    env_config = EnvConfig(
        Ssp(ssp_knots=ssp_knots, ssp_coefs=ssp_coefs, ssp_degree=ssp_degree),
        swell_height=0.0,
        bodies=bodies,
        halfspaces=halfspaces,
    )

    prog_config = ProgConfig(
        max_it=int(1e5), depth_step=1.0, save_to_csv=False, max_range=2e5
    )

    config = Config(prog_config=prog_config, env_config=env_config, sources=sources)

    config.write_to_json(save_path=f"{CONFIG_PATH}/run_config.json")
