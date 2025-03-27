from guacs import (
    # trace_rays,
    trace_beams,
    Config,
    ProgConfig,
    EnvConfig,
    SourceConfig,
    Body,
    # HalfSpace,
    Ssp,
)
import numpy as np
from scipy.interpolate import splev, splrep
import matplotlib.pyplot as plt


def munk_profile(depth: float) -> float:
    return 1500.0 * (
        1.0
        + 0.00737 * (2 * (depth - 1300) / 1500 - 1 + np.exp(-2 * (depth - 1300) / 1500))
    )


def fix_rays(raw_rays: list) -> list:
    for i in range(len(raw_rays)):
        raw_rays[i].range_vals = np.array(raw_rays[i].range_vals)
        raw_rays[i].depth_vals = np.array(raw_rays[i].depth_vals)
        raw_rays[i].time_vals = np.array(raw_rays[i].time_vals)
    return raw_rays


def plot_rays(cfg_file: Config, rays: list):
    bodies = cfg_file.env_config.bodies
    ssp = cfg_file.env_config.ssp
    knots = ssp.ssp_knots
    coefs = ssp.ssp_coefs
    degree = ssp.ssp_degree
    tck = (knots, coefs, degree)
    depth_range = list(range(0, 5000))
    c_profile = [splev(depth, tck, der=0) for depth in depth_range]
    g_profile = [splev(depth, tck, der=1) for depth in depth_range]

    fig, ax = plt.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [1, 4]})
    ax_ssp = ax[0]
    ax_ssp.invert_yaxis()
    ax_ssp.set_ylabel("Depth [$m$]")
    ax_ssp.set_xlabel("Sound speed [$m.s^{-1}$]")
    ax_ssp.plot(c_profile, depth_range, c="blue")
    ax_g = ax_ssp.twiny()
    ax_rays = ax[1]
    ax_g.set_xlabel("Sound gradient [$s^{-1}$]")
    ax_rays.set_xlabel("Range [$km$]")

    ax_g.plot(g_profile, depth_range, c="red")
    for bd in bodies:
        ax_rays.fill([rv / 1000 for rv in bd.range_vals], bd.depth_vals, c="r")

    for ray in rays:
        ax_rays.plot(
            np.array(ray.range_vals) / 1000, ray.depth_vals, c="k", linewidth=0.5
        )

    plt.show()


if __name__ == "__main__":
    sources = [
        SourceConfig(
            range_pos=0.0,
            depth_pos=1000.0,
            ray_fan_limits=(-0.2, 0.2),
            n_rays=100,
            source_level=150,
        )
    ]

    bodies = [
        Body(
            range_vals=[20000.0, 20100.0, 20100.0, 20000.0, 20000.0],
            depth_vals=[1400.0, 1400.0, 1310.0, 1310.0, 1400.0],
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

    # halfspaces = [
    #     HalfSpace(
    #         body=Body(
    #             range_vals=[0.0, 100000.0, 100000.0, 0.0, 0.0],
    #             depth_vals=[5000.0, 5000.0, 5100.0, 5100.0, 5000.0],
    #         ),
    #         sound_speed=1600.0,
    #         density=1.5,
    #     )
    # ]

    max_depth = 5000
    depth_step = 100
    munk_depths = range(-2 * depth_step, max_depth + 3 * depth_step, depth_step)
    munk_vals = [munk_profile(z) for z in munk_depths]

    (ssp_knots, ssp_coefs, ssp_degree) = splrep(munk_depths, munk_vals, k=3)

    env_config = EnvConfig(
        bodies=bodies,
        ssp=Ssp(ssp_knots=ssp_knots, ssp_coefs=ssp_coefs, ssp_degree=ssp_degree),
        swell_height=0.0,
        # halfspaces=halfspaces,
    )

    prog_config = ProgConfig(
        max_it=int(1e5),
        depth_step=1.0,
        max_range=2e5,
        min_range=-10.0,
        output_path="output_data",
    )

    config = Config(prog_config=prog_config, env_config=env_config, sources=sources)
    # rays = trace_rays(config)
    beams = trace_beams(config)
    rays = [bm.central_ray for bm in beams]
    rays = fix_rays(rays)
    plot_rays(config, rays)
