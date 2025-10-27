from guacs.guacs import (
    BeamResult,
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
import numpy as np

# import matplotlib.pyplot as plt
from typing import Tuple
import matplotlib
from python_utils import (
    munk_profile,
    # fix_rays,
    # plot_rays,
    # animate_propagation,
    # plot_pq,
    plot_sound_field,
)

matplotlib.use("QtAgg")

if __name__ == "__main__":
    sources = [
        SourceConfig(
            range_pos=0.0,
            depth_pos=1000.0,
            ray_fan_limits=(-0.1, 0.2),
            n_rays=10,
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

    # isospaces = [
    #     IsoSpace(
    #         body=Body(
    #             range_vals=[-100.0, 100000.0, 100000.0, -100.0, -100.0],
    #             depth_vals=[0.0, 0.0, 5000.0, 5000.0, 0.0],
    #         ),
    #         sound_speed=1600.0,
    #         density=1.5,
    #     )
    # ]

    isospaces = []

    max_depth = 5000
    depth_step = 1
    ssp_depths = range(-2 * depth_step, max_depth + 3 * depth_step, depth_step)
    ssp_vals = [munk_profile(z) for z in ssp_depths]

    spline_res: Tuple = splrep(ssp_depths, ssp_vals, k=3)

    env_config = EnvConfig(
        bodies=bodies,
        ssp=Ssp(
            ssp_knots=spline_res[0], ssp_coefs=spline_res[1], ssp_degree=spline_res[2]
        ),
        swell_height=0.0,
        isospaces=isospaces,
    )

    prog_config = ProgConfig(
        max_it=int(1e5),
        depth_step=1.0,
        max_range=5e4,
        min_range=-10.0,
        output_path="output_data",
    )

    ray_config = RayConfig(
        prog_config=prog_config, env_config=env_config, sources=sources
    )

    ranges = np.linspace(0, 5e4, 5001)
    depths = np.linspace(0, 5e3, 501)
    ranges, depths = np.meshgrid(ranges, depths)
    ranges = ranges.flatten()
    depths = depths.flatten()
    locs = [(ranges[i], depths[i]) for i in range(len(ranges))]

    beam_config = BeamConfig(
        ray_config=ray_config,
        pq_solver="Radau3IIA",
        pressure_locs=locs,
        window_width=None,
    )
    # rays = trace_rays(ray_config)
    beam_result: BeamResult = trace_beams(beam_config)
    plot_sound_field(cfg=beam_config.ray_config, beam_res=beam_result)
    # beams = beam_result.beams
    # rays = [bm.central_ray for bm in beams]
    # rays = fix_rays(rays)
    # ray_fig = plot_rays(ray_config, rays)
    # pq_fig = plot_pq(beams)
    # ani = animate_propagation(config, rays, 10, fast_fwd=100)
    # plt.show()
