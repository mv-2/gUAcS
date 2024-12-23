import matplotlib.pyplot as plt
from glob import glob
from python_interface_tools import Ray, Config
from typing import List
from scipy.interpolate import splev
import json


def plot_rays(rays: List[Ray]):
    cfg_file = json.load(open("configs/run_config.json"))
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

    for ray in rays:
        ax_rays.plot(ray.range_vals / 1000, ray.depth_vals, c="k")

    plt.show()


if __name__ == "__main__":
    ray_names = glob("output_data/rays/*.csv")
    rays = [None] * len(ray_names)
    for i, fn in enumerate(ray_names):
        rays[i] = Ray(fn)
    plot_rays(rays)
