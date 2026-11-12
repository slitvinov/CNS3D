import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ana

BLUE, ORANGE, AQUA, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#222222", "#888888"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6})


def bend(coarse, fine, out):
    P1 = ana.profiles(ana.read(coarse, ".o")[0])
    P8 = ana.profiles(ana.read(fine, ".o")[0])
    fig, ax = plt.subplots(2, 3, figsize=(9, 5.4), sharex=True, sharey=True)
    for r, (pn, zz) in enumerate(ana.PLANES):
        for c, th in enumerate(ana.STATIONS):
            a = ax[r, c]
            q = np.loadtxt(f"{ana.R}/fig19{pn}_theta{th}_humphrey.csv", delimiter=",", skiprows=1)
            p = np.loadtxt(f"{ana.R}/fig19{pn}_theta{th}_present.csv", delimiter=",", skiprows=1)
            o = np.argsort(p[:, 1])
            a.plot(p[o, 1], p[o, 0], color=ORANGE, lw=1.6, ls="--", label="report")
            a.plot(*P1[(th, pn)], color=AQUA, lw=1.6, label="this program, 20×17×9")
            a.plot(*P8[(th, pn)], color=BLUE, lw=2, label="this program, 153×129×65")
            a.plot(q[:, 1], q[:, 0], "o", ms=5, mfc="white", mec=INK, mew=1.2, label="Humphrey et al.")
            a.set_title(f"θ = {th}°, z/d = {zz:.2f}", fontsize=9, color=INK)
            a.set_xlim(0, 1)
            if r == 1:
                a.set_xlabel("(r − r_i)/d")
            if c == 0:
                a.set_ylabel("U/U_ref")
    ax[0, 0].legend(frameon=False, fontsize=8, loc="lower center")
    fig.tight_layout()
    fig.savefig(out, dpi=150)


def step(rows, out):
    fig, a = plt.subplots(figsize=(6, 4))
    for name, color, ds in rows:
        pts = sorted((1 / f, ana.xr(d)) for f, d in ds)
        x, y = zip(*pts)
        a.plot(x, y, "-o", color=color, lw=2, ms=6, mfc=color, mec="white", mew=1.5, label=name)
    for v, t in ((11.45, "Armaly et al. 11.45"), (10.53, "report, central 10.53"), (8.44, "report, hybrid 8.44")):
        a.axhline(v, color=MUTED, lw=1, ls=":")
        a.text(1.0, v, t, color=INK, fontsize=8, ha="right", va="bottom")
    a.set_xlim(0, 1.02)
    a.set_xlabel("x spacing / base x spacing")
    a.set_ylabel("reattachment length x_r / h")
    a.legend(frameon=False, fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(out, dpi=150)


if __name__ == "__main__":
    w, here = sys.argv[1], os.path.dirname(os.path.abspath(__file__))
    bend(f"{w}/bend1", f"{w}/bend8", f"{here}/bend.png")
    step((("y base", AQUA, ((1, f"{w}/it_lstep1_1e-6"), (4, f"{w}/ls_4_1"), (8, f"{w}/ls_8_1"), (16, f"{w}/ls_16_1"))),
          ("y refined 4 times", BLUE, ((4, f"{w}/ls_4_4"), (8, f"{w}/ls_8_4"), (16, f"{w}/ls_16_4"), (32, f"{w}/ls_32_4"), (40, f"{w}/ls_40_4"))),
          ("y refined 8 times", ORANGE, ((16, f"{w}/ls_16_8"), (32, f"{w}/ls_32_8")))), f"{here}/step.png")
