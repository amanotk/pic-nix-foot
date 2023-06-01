#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pathlib

import numpy as np
import matplotlib as mpl

mpl.use("Agg") if __name__ == "__main__" else None
from matplotlib import pyplot as plt

# global configuration
plt.rcParams.update({"font.size": 12})

if "PICNIX_DIR" in os.environ:
    sys.path.append(str(pathlib.Path(os.environ["PICNIX_DIR"]) / "script"))
import analysis


class Run(analysis.Run):
    def __init__(self, profile):
        super().__init__(profile)

    def summary2d(self, step):
        ## normalization
        parameter = self.config["parameter"]
        cc = parameter["cc"]
        wp = parameter["wp"]
        T = 1 / wp
        L = cc / wp

        data = self.read_field_at(step)
        uf = data["uf"]
        tt = self.get_field_time_at(step) / T
        xc = self.xc / L
        yc = self.yc / L
        xlim = (0, self.Nx * self.delh / L)
        ylim = (0, self.Ny * self.delh / L)

        ## figure and axes
        fig = plt.figure(1, figsize=(9.6, 7.2), dpi=120)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.08,
            left=0.10,
            right=0.95,
            hspace=0.10,
            wspace=0.35,
        )
        gridspec = fig.add_gridspec(5, 3, height_ratios=[1, 50, 20, 1, 50], width_ratios=[1, 1, 1])
        axs = [0] * 6
        cxs = [0] * 6
        for i in range(3):
            # for Ex, Ey, Ez
            axs[i + 0] = fig.add_subplot(gridspec[1, i])
            cxs[i + 0] = fig.add_subplot(gridspec[0, i])
            # for Bx, By, Bz
            axs[i + 3] = fig.add_subplot(gridspec[4, i])
            cxs[i + 3] = fig.add_subplot(gridspec[3, i])

        # coordinate
        X, Y = np.broadcast_arrays(xc[None, :], yc[:, None])

        # plot
        Ex = uf[..., 0].mean(axis=(0,))
        Ey = uf[..., 1].mean(axis=(0,))
        Ez = uf[..., 2].mean(axis=(0,))
        Bx = uf[..., 3].mean(axis=(0,))
        By = uf[..., 4].mean(axis=(0,))
        Bz = uf[..., 5].mean(axis=(0,))
        data = [Ex, Ey, Ez, Bx, By, Bz]
        name = [r"$E_x$", r"$E_y$", r"$E_z$", r"$B_x$", r"$B_y$", r"$B_z$"]

        for i in range(6):
            plt.sca(axs[i])
            axs[i].set_aspect("equal")
            # data
            plt.pcolormesh(X, Y, data[i], shading="nearest")
            # colorbar
            plt.colorbar(cax=cxs[i], orientation="horizontal")
            cxs[i].xaxis.set_ticks_position("top")
            cxs[i].set_title(name[i])
            # appearance
            axs[i].xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
            axs[i].yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
            axs[i].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))
            axs[i].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))
            axs[i].set_xlim(xlim)
            axs[i].set_ylim(ylim)
            axs[i].set_xlabel(r"$x / c/\omega_{pe}$")
            axs[i].set_ylabel(r"$y / c/\omega_{pe}$")
            ax_pos = axs[i].get_position()
            cx_pos = cxs[i].get_position()
            cxs[i].set_position([ax_pos.x0, cx_pos.y0, ax_pos.width, cx_pos.height])

        # title
        fig.suptitle(r"$\omega_{{pe}} t = {:6.2f}$".format(tt), x=0.5, y=0.99)

        return fig


def doit_job(profile, prefix, fps, cleanup):
    run = Run(profile)

    # for all snapshots
    for step in run.step_field:
        fig = run.summary2d(step)
        fig.savefig("{:s}-{:08d}.png".format(prefix, step))
        plt.close(fig)

    # convert to mp4
    analysis.convert_to_mp4("{:s}".format(prefix), fps, cleanup)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Summary Plot Script")
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="summary",
        help="Prefix used for output image and movie files",
    )
    parser.add_argument(
        "-f",
        "--fps",
        type=int,
        default=10,
        help="Frame/sec used for encoding movie file",
    )
    parser.add_argument(
        "-c",
        "--cleanup",
        action="store_true",
        default=False,
        help="Cleanup intermediate image files",
    )
    parser.add_argument("profile", nargs=1, help="run profile")

    args = parser.parse_args()
    profile = args.profile[0]
    doit_job(profile, args.prefix, args.fps, args.cleanup)
