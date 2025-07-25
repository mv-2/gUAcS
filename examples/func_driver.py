from guacs.guacs import (
    trace_rays,
    trace_beams,
    RayConfig,
    BeamConfig,
    ProgConfig,
    EnvConfig,
    SourceConfig,
    Body,
    IsoSpace,
    Ssp,
)
from scipy.interpolate import splrep
import matplotlib.pyplot as plt
from typing import Tuple
import matplotlib
from python_utils import munk_profile, fix_rays, plot_rays, animate_propagation, plot_pq

matplotlib.use("QtAgg")

if __name__ == "__main__":
    sources = [
        SourceConfig(
            range_pos=0.0,
            depth_pos=1000.0,
            ray_fan_limits=(-0.1, 0.2),
            n_rays=200,
            source_level=150,
            frequency=1000,
        )
    ]

    bodies = [
        # Body(
        #     range_vals=[20000.0, 20100.0, 20100.0, 20000.0, 20000.0],
        #     depth_vals=[1400.0, 1400.0, 1310.0, 1310.0, 1400.0],
        # ),
        Body(
            range_vals=[0.0, 200000.0, 200000.0, 0.0, 0.0],
            depth_vals=[5000.0, 5000.0, 5100.0, 5100.0, 5000.0],
        ),
        Body(
            range_vals=[0.0, 200000.0, 200000.0, 0.0, 0.0],
            depth_vals=[0.0, 0.0, -1.0, -1.0, 0.0],
        ),
    ]

    isospaces = [
        IsoSpace(
            body=Body(
                range_vals=[-100.0, 100000.0, 100000.0, -100.0, -100.0],
                depth_vals=[0.0, 0.0, 5000.0, 5000.0, 0.0],
            ),
            sound_speed=1600.0,
            density=1.5,
        )
    ]

    max_depth = 5000
    depth_step = 100
    ssp_depths = range(-2 * depth_step, max_depth + 3 * depth_step, depth_step)
    ssp_vals = [munk_profile(z) for z in ssp_depths]

    spline_res: Tuple = splrep(ssp_depths, ssp_vals, k=3)

    env_config = EnvConfig(
        bodies=bodies,
        ssp=Ssp(
            ssp_knots=spline_res[0], ssp_coefs=spline_res[1], ssp_degree=spline_res[2]
        ),
        swell_height=0.0,
        isospaces=[],
    )

    prog_config = ProgConfig(
        max_it=int(100),
        depth_step=1.0,
        max_range=2e5,
        min_range=-10.0,
        output_path="output_data",
    )

    ray_config = RayConfig(
        prog_config=prog_config, env_config=env_config, sources=sources
    )
    beam_config = BeamConfig(
        ray_config=ray_config,
        pq_solver="Radau3IIA",
        pressure_locs=[],
        window_width=1000,
    )
    rays = trace_rays(ray_config)
    # beams = trace_beams(beam_config)
    # rays = [bm.central_ray for bm in beams]
    rays = fix_rays(rays)
    ray_fig = plot_rays(ray_config, rays)
    # pq_fig = plot_pq(beams)
    # ani = animate_propagation(config, rays, 10, fast_fwd=100)
    plt.show()
