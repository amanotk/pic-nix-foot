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

# result of linear analysis
LINEAR_THEORY_NPZ = "alfven-linear-theory.npz"


class Run(analysis.Run):
    def __init__(self, profile):
        super().__init__(profile)

        # read and store all data
        Nstep = len(self.file_field)
        parameter = self.config["parameter"]
        Nx = parameter["Nx"]
        Ny = parameter["Ny"]
        Nz = parameter["Nz"]
        delt = parameter["delt"]
        delh = parameter["delh"]
        cc = parameter["cc"]
        wp = parameter["wp"]
        mime = parameter["mime"]
        sigma = parameter["sigma"]
        vai = cc * np.sqrt(sigma / mime)
        wci = wp * np.sqrt(sigma) / mime
        tt = np.zeros((Nstep,), dtype=np.float64)
        uf = np.zeros((Nstep, Nz, Ny, Nx, 6), dtype=np.float64)
        um = np.zeros((Nstep, Nz, Ny, Nx, 2, 11), dtype=np.float64)

        # read
        for i in range(Nstep):
            step = self.step_field[i]
            data = self.read_field_at(step)
            tt[i] = self.get_field_time_at(step)
            uf[i, ...] = data["uf"]
            um[i, ...] = data["um"]

        # store as dict
        self.data = dict(
            Nx=Nx,
            Ny=Ny,
            Nz=Nz,
            delt=delt,
            delh=delh,
            cc=cc,
            wp=wp,
            mime=mime,
            sigma=sigma,
            vai=vai,
            wci=wci,
            xc=self.xc,
            tt=tt,
            uf=uf,
            um=um,
        )

    def linear_theory(self):
        Nx = self.data["Nx"]
        delh = self.data["delh"]
        vai = self.data["vai"]
        wci = self.data["wci"]
        tt = self.data["tt"] * wci
        uf = self.data["uf"]

        # take Fourier transform for transverse B-field in complex representation
        bt = (uf[..., 4] - 1j * uf[..., 5]).mean(axis=(1, 2))
        kk = np.fft.fftfreq(Nx, delh / (vai / wci)) * 2 * np.pi
        bk = np.fft.fft(bt, axis=-1) / Nx

        # result of linear dispersion analysis
        disp = np.load(LINEAR_THEORY_NPZ)

        tzero = 10.0
        mode_to_plot = [2, 3, 4, 5, 6]

        # plot
        fig, axs = plt.subplots(len(mode_to_plot), figsize=(6, 10), sharex=True)
        fig.subplots_adjust(
            top=0.95,
            bottom=0.05,
            left=0.15,
            right=0.95,
            hspace=0.20,
            wspace=0.10,
        )

        for imode, mode in enumerate(mode_to_plot):
            plt.sca(axs[imode])
            # plot simulation result
            ik = disp["k"].searchsorted(kk[mode])
            gm = disp["wim"][ik]
            plt.plot(tt, np.abs(bk[:, +mode]), "b-", label="simulation")
            # plot theoretical growth curve
            izero = tt.searchsorted(tzero)
            bzero = np.abs(bk[izero, +mode])
            plt.plot(tt, bzero * np.exp(gm * (tt - tzero)), "r--", label="theory")
            # appearance
            plt.title(r"mode = {:d} $(k v_{{A}}/\Omega_{{ci}}$ = {:5.2f})".format(mode, kk[mode]))
            plt.semilogy()
            plt.legend(loc="lower right")
            plt.ylabel(r"$|B_y - i B_z|/B_0$".format(mode))
            plt.ylim(1e-4, 1e0)
            plt.grid()

        plt.sca(axs[-1])
        plt.xlabel(r"$\Omega_{ci} t$")

        return fig

    def energy_history(self):
        wci = self.data["wci"]
        cc = self.data["cc"]
        tt = self.data["tt"] * wci
        uf = self.data["uf"]
        um = self.data["um"]

        # show energy history
        Wf = 0.5 * np.sum(uf**2, axis=(1, 2, 3, 4))
        We = np.sum(um[..., 0, 4] * cc - um[..., 0, 0] * cc**2, axis=(1, 2, 3))
        Wi = np.sum(um[..., 1, 4] * cc - um[..., 1, 0] * cc**2, axis=(1, 2, 3))
        Wtot = Wf[0] + We[0] + Wi[0]

        fig, axs = plt.subplots(2, figsize=(8, 6), sharex=True)
        fig.subplots_adjust(
            top=0.95,
            bottom=0.15,
            left=0.15,
            right=0.80,
            hspace=0.10,
            wspace=0.10,
        )

        plt.sca(axs[0])
        plt.plot(tt, Wf / Wtot, label="field")
        plt.plot(tt, We / Wtot, label="electron")
        plt.plot(tt, Wi / Wtot, label="ion")
        plt.ylabel("Energy")
        plt.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        plt.grid()

        plt.sca(axs[1])
        plt.plot(tt, 100 * ((Wf + We + Wi) / Wtot - 1), "k-")
        plt.xlabel(r"$\Omega_{ci} t$")
        plt.ylabel("Energy Conservation Error [%]")
        plt.xlim(tt[0], tt[-1])
        plt.grid()

        for ax in axs:
            ax.set_axisbelow(True)

        return fig

    def helicity_decomposition(self):
        Nx = self.data["Nx"]
        vai = self.data["vai"]
        wci = self.data["wci"]
        xc = self.data["xc"] / (vai / wci)
        tt = self.data["tt"] * wci
        uf = self.data["uf"]

        # show time evolution with helicity decomposition
        by = uf[..., 4].mean(axis=(1, 2))
        bz = uf[..., 5].mean(axis=(1, 2))
        bt = by - 1j * bz
        bk = np.fft.fft(bt, axis=-1)
        bp = np.zeros_like(bk)
        bm = np.zeros_like(bk)

        ip = np.arange(0, Nx // 2 + 1, 1)
        im = np.arange(Nx, Nx // 2 - 1, -1)
        im[0] = 0

        # positive helicity
        bp[...] = 0.0
        bp[..., ip] = bk[..., ip]
        bp = np.fft.ifft(bp, axis=-1)
        bpy = +bp.real
        bpz = -bp.imag

        # positive helicity
        bm[...] = 0.0
        bm[..., im] = bk[..., im]
        bm = np.fft.ifft(bm, axis=-1)
        bmy = +bm.real
        bmz = -bm.imag

        # plot
        fig = plt.figure(1, figsize=(8, 4), dpi=120)
        fig.subplots_adjust(
            top=0.85,
            bottom=0.15,
            left=0.10,
            right=0.95,
            hspace=0.10,
            wspace=0.20,
        )
        gridspec = fig.add_gridspec(
            2, 3, height_ratios=[2, 50], width_ratios=[1, 1, 1], hspace=0.05, wspace=0.25
        )

        axs = [0] * 3
        cxs = [0] * 3
        for i in range(3):
            axs[i] = fig.add_subplot(gridspec[1, i])
            cxs[i] = fig.add_subplot(gridspec[0, i])

        X, T = np.broadcast_arrays(xc[None, :], tt[:, None])

        plt.sca(axs[0])
        plt.pcolormesh(X, T, bz, shading="nearest", clim=[-0.8, +0.8])
        plt.colorbar(cax=cxs[0], orientation="horizontal")
        cxs[0].xaxis.set_ticks_position("top")
        cxs[0].set_title(r"$B_z$ (raw)")

        plt.sca(axs[1])
        plt.pcolormesh(X, T, bpz, shading="nearest", clim=[-0.8, +0.8])
        plt.colorbar(cax=cxs[1], orientation="horizontal")
        cxs[1].xaxis.set_ticks_position("top")
        cxs[1].set_title(r"$B_z^{+}$")

        plt.sca(axs[2])
        plt.pcolormesh(X, T, bmz, shading="nearest", clim=[-0.1, +0.1])
        plt.colorbar(cax=cxs[2], orientation="horizontal")
        cxs[2].xaxis.set_ticks_position("top")
        cxs[2].set_title(r"$B_z^{-}$")

        for ax in axs:
            ax.set_xlabel(r"$x / c/\omega_{pi}$")
            ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1.0))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1.0))
        axs[0].set_ylabel(r"$\Omega_{ci} t$")

        return fig

    def velocity_distribution(self, step, vmin=None, vmax=None):
        vai = self.data["vai"]
        wci = self.data["wci"]

        # read and plot velocity distribution function
        up = self.read_particle_at(step)
        tp = self.get_particle_time_at(step) * wci

        vx = up[1][:, 3] / vai
        vy = up[1][:, 4] / vai
        hist2d = analysis.Histogram2D(vx, vy, (-2.0, +4.0, 151), (-3.0, +3.0, 151))
        X, Y, Z = hist2d.pcolormesh_args()
        vmax = 1.0e-2 if vmax is None else vmax
        vmin = vmax * 1.0e-3 if vmin is None else vmin
        norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

        # plot
        cb_ratio = 0.05
        figx = 4 * (1 + 9 * cb_ratio)
        figy = 4 * (1 + 4 * cb_ratio)
        fig = plt.figure(1, figsize=(figx, figy), dpi=120)
        fig.subplots_adjust(
            top=1 - 2 * cb_ratio,
            bottom=2 * cb_ratio,
            left=3 * cb_ratio,
            right=1 - 4 * cb_ratio,
            hspace=cb_ratio,
            wspace=cb_ratio,
        )
        gridspec = fig.add_gridspec(
            1,
            2,
            height_ratios=[1],
            width_ratios=[1, cb_ratio],
        )

        ax = fig.add_subplot(gridspec[0, 0])
        cx = fig.add_subplot(gridspec[0, 1])

        plt.sca(ax)
        plt.pcolormesh(X, Y, Z, norm=norm)
        ax = plt.gca()
        ax.set_aspect("equal")
        ax.set_xlim(-2.0, +4.0)
        ax.set_ylim(-3.0, +3.0)
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
        ax.grid()
        ax.set_axisbelow(True)
        ax.set_xlabel(r"$v_x / v_{A}$")
        ax.set_ylabel(r"$v_y / v_{A}$")
        ax.set_title(r"$\Omega_{{ci}} t =$ {:5.2f}".format(tp))

        # colorbar
        plt.colorbar(cax=cx, orientation="vertical")
        cx.set_ylabel(r"$f(v)$")

        ax_pos = ax.get_position()
        cx_pos = cx.get_position()
        cx.set_position([cx_pos.x1, ax_pos.y0, cx_pos.width, ax_pos.height])

        return fig


def doit_job(profile, prefix, fps, cleanup):
    run = Run(profile)

    # comparison with linear theory
    fig = run.linear_theory()
    fig.savefig("{:s}-linear.png".format(prefix))
    plt.close(fig)

    # energy history
    fig = run.energy_history()
    fig.savefig("{:s}-energy.png".format(prefix))
    plt.close(fig)

    # time evolution with helicity decomposition
    fig = run.helicity_decomposition()
    fig.savefig("{:s}-helicity.png".format(prefix))
    plt.close(fig)

    # for all snapshots
    for step in run.step_particle:
        fig = run.velocity_distribution(step)
        fig.savefig("{:s}-vdf-{:08d}.png".format(prefix, step))
        plt.close(fig)
    # convert to mp4
    analysis.convert_to_mp4("{:s}-vdf".format(prefix), fps, cleanup)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Quicklook Script")
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="alfven",
        help="Prefix used for output image and movie files",
    )
    parser.add_argument(
        "-f",
        "--fps",
        type=int,
        default=2,
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
