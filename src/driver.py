from guacs import run_sim
import numpy as np
from python_interface_tools import Body
from scipy.interpolate import splev
import json
import matplotlib.pyplot as plt

CONFIG_PATH = "configs/run_config.json"
OUTPUT_PATH = "output_data"


def fix_rays(raw_rays: list) -> list:
    for i in range(len(raw_rays)):
        raw_rays[i].range_vals = np.array(raw_rays[i].range_vals)
        raw_rays[i].depth_vals = np.array(raw_rays[i].depth_vals)
        raw_rays[i].time_vals = np.array(raw_rays[i].time_vals)
    return raw_rays


def plot_rays(rays: list):
    cfg_file = json.load(open("configs/run_config.json"))
    bodies = []
    for bd in cfg_file["env_config"]["bodies"]:
        bodies += [Body(**bd)]
    ssp = cfg_file["env_config"]["ssp"]
    knots = ssp["ssp_knots"]
    coefs = ssp["ssp_coefs"]
    degree = ssp["ssp_degree"]
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
    rays = run_sim(CONFIG_PATH, OUTPUT_PATH)
    rays = fix_rays(rays)
    plot_rays(rays)
