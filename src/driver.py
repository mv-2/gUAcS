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
    # Ray,
    PyBeam,
)
import numpy as np
from scipy.interpolate import splev, splrep
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


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


def plot_environment(cfg_file: Config) -> list:
    fig, ax = plt.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [1, 4]})
    bodies = cfg_file.env_config.bodies
    ssp = cfg_file.env_config.ssp
    knots = ssp.ssp_knots
    coefs = ssp.ssp_coefs
    degree = ssp.ssp_degree
    tck = (knots, coefs, degree)
    depth_range = list(range(0, 5000))
    c_profile = [splev(depth, tck, der=0) for depth in depth_range]
    g_profile = [splev(depth, tck, der=1) for depth in depth_range]

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
    ax_rays.set_xlim(
        [cfg_file.prog_config.min_range / 1000, cfg_file.prog_config.max_range / 1000]
    )
    return [fig, ax_rays, ax_ssp, ax_g]


def plot_pq(beam: PyBeam):
    fig, ax = plt.subplots(2, 2)
    iters = list(range(len(beam.p_re)))
    ax[0, 0].plot(iters, beam.p_re)
    ax[0, 0].set_ylabel("Re($p$)")
    ax[0, 0].set_xlabel("Iterations")
    ax[0, 0].set_ylim([-1, 1e5])

    ax[1, 0].plot(iters, beam.p_im)
    ax[1, 0].set_ylabel("Im($p$)")
    ax[1, 0].set_xlabel("Iterations")
    ax[1, 0].set_ylim([-1, 1e5])

    ax[0, 1].plot(iters, beam.q_re)
    ax[0, 1].set_ylabel("Re($q$)")
    ax[0, 1].set_xlabel("Iterations")
    ax[0, 1].set_ylim([-1, 1e5])

    ax[1, 1].plot(iters, beam.q_im)
    ax[1, 1].set_ylabel("Im($q$)")
    ax[1, 1].set_xlabel("Iterations")
    ax[1, 1].set_ylim([-1, 1e5])

    plt.show()


def plot_rays(cfg_file: Config, rays: list):
    fig, ax_rays, ax_ssp, ax_g = plot_environment(cfg_file)

    for ray in rays:
        ax_rays.plot(
            np.array(ray.range_vals) / 1000, ray.depth_vals, c="k", linewidth=0.5
        )

    plt.show()


def time_interp_rays(rays: list, time_vals: np.array) -> list:
    for i in range(len(rays)):
        rays[i].range_vals = np.interp(
            time_vals, rays[i].time_vals, rays[i].range_vals, right=np.nan
        )
        rays[i].depth_vals = np.interp(
            time_vals, rays[i].time_vals, rays[i].depth_vals, right=np.nan
        )
        rays[i].time_vals = time_vals
    return rays


def animate_propagation(cfg_file: Config, rays: list, framerate: float):
    fig, ax_rays, ax_ssp, ax_g = plot_environment(cfg_file)
    max_time = np.nanmax([np.nanmax(ray.time_vals) for ray in rays])
    time_vals = np.arange(0.0, max_time, 1 / framerate)
    rays = time_interp_rays(rays, time_vals)
    n_rays = len(rays)
    lines = [None] * n_rays
    for i in range(n_rays):
        (lines[i],) = ax_rays.plot([], [], c="k", linewidth=0.5)
        rays[i].range_vals = np.array(rays[i].range_vals) / 1000

    def update_frame(i):
        for j in range(n_rays):
            lines[j].set_xdata(rays[j].range_vals[:i])
            lines[j].set_ydata(rays[j].depth_vals[:i])

    _ = FuncAnimation(
        fig, update_frame, frames=len(time_vals), interval=1000 / framerate
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
            frequency=500,
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
        depth_step=10.0,
        max_range=2e5,
        min_range=-10.0,
        output_path="output_data",
        pq_solver="Radau3IIA",
    )

    config = Config(prog_config=prog_config, env_config=env_config, sources=sources)
    # rays = trace_rays(config)
    beams = trace_beams(config)
    rays = [bm.central_ray for bm in beams]
    rays = fix_rays(rays)
    plot_rays(config, rays)
    plot_pq(beams[0])
    # animate_propagation(config, rays, 10)
