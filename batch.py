#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

os.environ["HDF5_EXTFILE_PREFIX"] = "${ORIGIN}"

import json
import tqdm
import glob
import concurrent
import concurrent.futures
import numpy as np
import h5py
import matplotlib as mpl

# for batch mode
mpl.use("Agg")
from matplotlib import pyplot as plt

summary_plot_config = dict(
    figure=dict(
        width=9.6,
        height=3.6,
    ),
    subplot=dict(top=0.80, bottom=0.14, left=0.075, right=0.95, hspace=0.05, wspace=0.35),
)

# command line for ffmpeg
ffmpeg_cmdline_format = (
    "ffmpeg -r {fps:d} -i {src:s} -vcodec libx264 -pix_fmt yuv420p -r 60 {dst:s}"
)


class Run(object):
    def __init__(self, cfgfile):
        self.dirname = os.path.dirname(cfgfile)
        self.read_config(cfgfile)
        self.read_coord(self.file_field)
        self.read_chunkmap()

    def read_config(self, cfgfile):
        cfg = json.loads(open(cfgfile, "r").read())
        self.cfg = cfg
        self.Ns = cfg["parameter"]["Ns"]
        self.Nx = cfg["parameter"]["Nx"]
        self.Ny = cfg["parameter"]["Ny"]
        self.Nz = cfg["parameter"]["Nz"]
        self.Cx = cfg["parameter"]["Cx"]
        self.Cy = cfg["parameter"]["Cy"]
        self.Cz = cfg["parameter"]["Cz"]
        self.delt = cfg["parameter"]["delt"]
        self.delh = cfg["parameter"]["delh"]
        for diagnostic in cfg["diagnostic"]:
            prefix = diagnostic["prefix"]
            path = os.sep.join([self.dirname, diagnostic["path"]])
            interval = diagnostic["interval"]
            file = sorted(glob.glob(os.sep.join([path, prefix]) + "*.h5"))
            # field
            if diagnostic["name"] == "field":
                self.file_field = file
                self.time_field = np.arange(len(file)) * self.delt * interval
                self.data_field = {}
            # particle
            if diagnostic["name"] == "particle":
                self.file_particle = file
                self.time_particle = np.arange(len(file)) * self.delt * interval
                self.data_particle = {}

    def read_coord(self, files):
        with h5py.File(files[0], "r") as h5fp:
            self.xc = np.unique(h5fp.get("xc")[()])
            self.yc = np.unique(h5fp.get("yc")[()])
            self.zc = np.unique(h5fp.get("zc")[()])

    def read_chunkmap(self):
        step = 0
        with h5py.File(self.file_field[step], "r") as h5fp:
            chunkid = h5fp.get("/chunkmap/chunkid")[()]
            rank = h5fp.get("/chunkmap/rank")[()]
            coord = h5fp.get("/chunkmap/coord")[()]
        cdelx = self.Nx // self.Cx * self.delh
        cdely = self.Ny // self.Cy * self.delh
        cdelz = self.Nz // self.Cz * self.delh
        self.chunkmap = {
            "chunkid": chunkid,
            "rank": rank,
            "coord": coord,
            "cdelx": cdelx,
            "cdely": cdely,
            "cdelz": cdelz,
        }

    def read_field_at(self, step):
        if self.data_field.get("step", -1) == step:
            # return cache
            return [self.data_field.get(var) for var in ("uf", "uj", "um")]
        else:
            Ns = self.Ns
            Nz = self.Nz
            Ny = self.Ny
            Nx = self.Nx
            with h5py.File(self.file_field[step], "r") as h5fp:
                uf = h5fp.get("/vds/uf")[()]
                uj = h5fp.get("/vds/uf")[()]
                um = h5fp.get("/vds/um")[()]
            # cache
            self.data_field = {
                "step": step,
                "uf": uf,
                "uj": uj,
                "um": um,
            }
        return uf, uj, um


def convert_to_mp4(prefix, fps, cleanup):
    import os
    import glob
    import subprocess

    src = prefix + "-%08d.png"
    dst = prefix + ".mp4"
    cmd = ffmpeg_cmdline_format.format(src=src, dst=dst, fps=fps)

    # run
    status = True
    try:
        if os.path.exists(dst):
            os.remove(dst)
        subprocess.run(cmd.split(), capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        status = False
        print(e.cmd)
        print(e.returncode)
        print(e.output)
        print(e.stdout)
        print(e.stderr)

    # cleanup if succeeded
    if status and cleanup:
        # remove files
        files = list()
        files.extend(glob.glob("{:s}-*.png".format(prefix)))
        for f in files:
            os.remove(f)


def doit_job(cfgfile, step, prefix):
    run = Run(cfgfile)
    output = "{:s}-{:08d}.png".format(prefix, step)
    uf, uj, um = run.read_field_at(step)
    summary_plot(run, uf, step, output)
    return True


def doit_parallel(cfgfile, prefix, numprocess, cleanup):
    run = Run(cfgfile)
    numstep = len(run.file_field)

    # parallel execution with progress bar
    with concurrent.futures.ProcessPoolExecutor(max_workers=numprocess) as executor:
        future_list = []
        for step in range(numstep):
            f = executor.submit(doit_job, cfgfile, step, prefix)
            future_list.append(f)
        # show progress
        progress_bar = tqdm.tqdm(total=numstep, desc="Generating Plot")
        progress_bar.update(0)
        for f in concurrent.futures.as_completed(future_list):
            progress_bar.update(1)

    # convert to mp4
    convert_to_mp4(prefix, 10, cleanup)


def summary_plot(run, uf, step, output=None):
    config = summary_plot_config

    ## ion-scale normalizeation
    parameter = run.cfg["parameter"]
    cc = parameter["cc"]
    wp = parameter["wp"]
    mime = parameter["mime"]
    T = 1 / (wp / np.sqrt(mime))
    L = cc / wp * np.sqrt(mime)

    xlim = [0, run.Nx * run.delh / L]
    ylim = [0, run.Ny * run.delh / L]
    time = run.time_field[step] / T
    xc = run.xc / L
    yc = run.yc / L

    ## figure and axes
    figw = config["figure"]["width"]
    figh = config["figure"]["height"]
    figure = plt.figure(1, figsize=(figw, figh), dpi=120)
    figure.clf()
    gridspec = figure.add_gridspec(
        2, 3, width_ratios=[1, 1, 1], height_ratios=[2, 50], **config["subplot"]
    )
    axs = list()
    for i in range(3):
        ax = figure.add_subplot(gridspec[1, i])
        ax.set_aspect("equal")
        axs.append(ax)

    # coordinate
    X, Y = np.broadcast_arrays(xc[None, :], yc[:, None])

    # Bx, By, Bz
    Bx = uf[..., 3].mean(axis=(0,))
    By = uf[..., 4].mean(axis=(0,))
    Bz = uf[..., 5].mean(axis=(0,))

    ## Bx
    plt.sca(axs[0])
    plt.pcolormesh(X, Y, Bx, shading="nearest")
    cax = figure.add_subplot(gridspec[0, 0])
    plt.colorbar(cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.set_title(r"$B_x$")

    ## By
    plt.sca(axs[1])
    plt.pcolormesh(X, Y, By, shading="nearest")
    cax = figure.add_subplot(gridspec[0, 1])
    plt.colorbar(cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.set_title(r"$B_y$")

    ## Bz
    plt.sca(axs[2])
    plt.pcolormesh(X, Y, Bz, shading="nearest")
    cax = figure.add_subplot(gridspec[0, 2])
    plt.colorbar(cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("top")
    cax.set_title(r"$B_z$")

    ## appearance
    for ax in axs:
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(5.0))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1.0))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(5.0))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1.0))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(r"$x / c/\omega_{pi}$")
        ax.set_ylabel(r"$y / c/\omega_{pi}$")
    figure.suptitle(r"$\omega_{pi} t$" + " = {:6.2f}".format(time), x=0.5, y=0.98)

    ## save
    if output is not None:
        plt.savefig(output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Quicklook Script")
    parser.add_argument(
        "-n",
        "--numprocess",
        type=int,
        default=4,
        help="Number of processes to use for parallel execution",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="foot",
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
    parser.add_argument("config", nargs=1, help="configuration file")

    args = parser.parse_args()
    config = args.config[0]
    prefix = args.prefix
    numprocess = args.numprocess
    cleanup = args.cleanup
    doit_parallel(config, prefix, numprocess, cleanup)
