from guacs.guacs import RayConfig, Ray, PyBeam, BeamConfig, BeamResult
import numpy as np
from scipy.interpolate import splev
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure
from matplotlib.artist import Artist
from typing import List, Iterable


def munk_profile(depth: float) -> float:
    # Munk sound speed profile defined by depth in m
    return 1500.0 * (
        1.0
        + 0.00737 * (2 * (depth - 1300) / 1500 - 1 + np.exp(-2 * (depth - 1300) / 1500))
    )


def downward_refract(depth: float) -> float:
    # Downward refraction sound speed profile defined by depth in m
    return 1677.3319 / np.sqrt(1 + 2.0 * 1.2276762 * depth / 1677.3319)


def fix_rays(raw_rays: list) -> List[Ray]:
    # Sets all ray values to array rather than list
    for i in range(len(raw_rays)):
        raw_rays[i].range_vals = np.array(raw_rays[i].range_vals)
        raw_rays[i].depth_vals = np.array(raw_rays[i].depth_vals)
        raw_rays[i].time_vals = np.array(raw_rays[i].time_vals)
    return raw_rays


def plot_environment(cfg_file: RayConfig) -> List:
    # Plots 2D simulation environment based on RayConfig object. Returns figure and axes handles
    fig, ax = plt.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [1, 4]})
    bodies = cfg_file.env_config.bodies
    ssp = cfg_file.env_config.ssp
    isospaces = cfg_file.env_config.isospaces
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
    ax_g.plot(g_profile, depth_range, c="red")
    ax_g.set_xlabel("Sound gradient [$s^{-1}$]")

    ax_rays = ax[1]
    for bd in bodies:
        ax_rays.fill([rv / 1000 for rv in bd.range_vals], bd.depth_vals, c="r")

    for ispc in isospaces:
        ax_rays.fill(
            [rv / 1000 for rv in ispc.body.range_vals], ispc.body.depth_vals, c="orange"
        )

    ax_rays.set_xlabel("Range [$km$]")
    ax_rays.set_xlim(
        [cfg_file.prog_config.min_range / 1000, cfg_file.prog_config.max_range / 1000]
    )

    return [fig, ax_rays, ax_ssp, ax_g]


def plot_pq(beams: List[PyBeam]) -> Figure:
    # plot evolution of pq derivatives
    fig, ax = plt.subplots(2, 2)
    for beam in beams:
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

    return fig


def plot_rays(cfg: RayConfig, rays: list) -> Figure:
    # plot resultant rays on backdrop defined by RayConfig

    fig, ax_rays, _, _ = plot_environment(cfg)

    for ray in rays:
        ax_rays.plot(
            np.array(ray.range_vals) / 1000, ray.depth_vals, c="k", linewidth=0.5
        )

    return fig


def time_interp_rays(rays: list, time_vals: np.ndarray) -> List[Ray]:
    # resample ray positions in equispaced time values required for animating propagation
    for i in range(len(rays)):
        rays[i].range_vals = np.interp(
            time_vals, rays[i].time_vals, rays[i].range_vals, right=np.nan
        )
        rays[i].depth_vals = np.interp(
            time_vals, rays[i].time_vals, rays[i].depth_vals, right=np.nan
        )
        rays[i].time_vals = time_vals
    return rays


def animate_propagation(
    cfg: RayConfig, rays: List[Ray], framerate: float, fast_fwd: float = 1
):
    # Creates ray propagation animation. Note that timesteps must be sufficiently small to ensure that ray follows smooth path.
    fig, ax_rays, _, _ = plot_environment(cfg)
    max_time = np.nanmax([np.nanmax(ray.time_vals) for ray in rays])
    time_vals = np.arange(0.0, max_time, 1 / framerate)
    rays = time_interp_rays(rays, time_vals)
    n_rays = len(rays)
    lines = np.zeros([n_rays], dtype=object)
    for i in range(n_rays):
        (lines[i],) = ax_rays.plot([], [], c="k", linewidth=0.5)
        rays[i].range_vals = np.array(rays[i].range_vals) / 1000

    def update_frame(i: int) -> Iterable[Artist]:
        for j in range(n_rays):
            lines[j].set_xdata(rays[j].range_vals[:i])
            lines[j].set_ydata(rays[j].depth_vals[:i])
        return lines

    ani = FuncAnimation(
        fig, update_frame, frames=len(time_vals), interval=1000 / (framerate * fast_fwd)
    )
    return ani


def plot_sound_field(cfg: BeamConfig, beam_res: BeamResult):
    # Plot heatmap of sound pressure
    _, ax, _, _ = plot_environment(cfg)
    ranges = [loc[0] / 1000.0 for loc in beam_res.pressures.locations]
    depths = [loc[1] for loc in beam_res.pressures.locations]
    magnitudes = beam_res.pressures.mag

    n_ranges = len(np.unique(ranges))
    n_depths = len(np.unique(depths))
    new_shape = (n_ranges, n_depths)
    ranges = np.reshape(ranges, new_shape)
    depths = np.reshape(depths, new_shape)
    magnitudes = np.reshape(magnitudes, new_shape)

    ax.pcolormesh(
        ranges,
        depths,
        magnitudes,
        shading="gouraud",
        cmap="jet",
    )
    # plt.colorbar(None, cax=ax)

    plt.show()
