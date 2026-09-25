#!/usr/bin/env python3
"""Build every figure of the report from results/*.txt and build/*.bin -> report/fig/*.png"""
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kimg  # noqa: E402

ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
BUILD = os.path.join(ROOT, "build")
OUT = os.path.join(ROOT, "report", "fig")
os.makedirs(OUT, exist_ok=True)

# reference categorical palette (fixed order), text/grid tokens
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.family": ["DejaVu Sans", "WenQuanYi Zen Hei"],
    "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "lines.linewidth": 2, "savefig.dpi": 160, "savefig.bbox": "tight",
})


def read(name):
    with open(os.path.join(RES, name)) as f:
        return f.read()


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- renders
def fig_render():
    d = kimg.load(os.path.join(BUILD, "render_double.bin"))
    s = kimg.load(os.path.join(BUILD, "render_shadow.bin"))
    ext = [-20, 20, -20, 20]
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.8))
    ax[0].imshow(d["I"] ** 0.5, origin="lower", cmap="afmhot", extent=ext)
    ax[0].set_title("thin disk, a=0.9375, i=60° (double, RK4)")
    ax[1].imshow(s["status"] == 1, origin="lower", cmap="Greys", extent=ext)
    # analytic Bardeen curve
    a, inc = 0.9375, np.radians(60)
    rp, rm = 2 * (1 + np.cos(2 / 3 * np.arccos(-a))), 2 * (1 + np.cos(2 / 3 * np.arccos(a)))
    r = np.linspace(rp, rm, 4000)
    xi = (r**2 * (3 - r) - a**2 * (r + 1)) / (a * (r - 1))
    eta = r**3 * (4 * a**2 - r * (r - 3) ** 2) / (a**2 * (r - 1) ** 2)
    b2 = eta + a**2 * np.cos(inc) ** 2 - xi**2 / np.tan(inc) ** 2
    ok = b2 >= 0
    al, be = -xi[ok] / np.sin(inc), np.sqrt(b2[ok])
    ax[1].plot(al, be, color=SERIES[1], lw=1.5, label="Bardeen critical curve")
    ax[1].plot(al, -be, color=SERIES[1], lw=1.5)
    ax[1].set_xlim(-8, 8), ax[1].set_ylim(-8, 8)
    ax[1].legend(loc="upper right", fontsize=7)
    ax[1].set_title("captured rays (black) vs analytic shadow")
    im = ax[2].imshow(d["steps"], origin="lower", cmap="Blues", extent=ext)
    fig.colorbar(im, ax=ax[2], fraction=0.046, label="RK4 steps per ray")
    ax[2].set_title("cost per ray (load imbalance)")
    for x in ax:
        x.set_xlabel("α [M]"), x.set_ylabel("β [M]"), x.grid(False)
    save(fig, "render.png")


def fig_raptor():
    def load(f):
        with open(f) as h:
            head = h.readline().split()
            v = np.array([float(x) for x in h.read().split()])
        n = int(head[0])
        return v.reshape(n, n).T  # file order: i (x) outer, j (y) inner
    a = load(os.path.join(BUILD, "raptor_O3/output/img_data_0_1.000000e+11_60.00.dat"))
    b = load(os.path.join(BUILD, "raptor_Ofast/output/img_data_0_1.000000e+11_60.00.dat"))
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.8))
    ax[0].imshow(a, origin="lower", cmap="afmhot", extent=[-20, 20, -20, 20])
    ax[0].set_title("RAPTOR GRMHD image (dump040, 100 GHz), -O3")
    rel = np.where(a > 0, np.abs(a - b) / np.maximum(a, 1e-300), 0)
    ys, xs = np.nonzero(rel)
    ax[1].imshow(np.zeros_like(a), origin="lower", cmap="Greys", extent=[-20, 20, -20, 20], vmin=0, vmax=1)
    sc = ax[1].scatter(-20 + (xs + 0.5) * 40 / a.shape[1], -20 + (ys + 0.5) * 40 / a.shape[0], c=np.log10(rel[ys, xs]),
                       s=14, cmap="Oranges", vmin=-6, vmax=-1, edgecolors="none")
    fig.colorbar(sc, ax=ax[1], fraction=0.046, label="log10 |Ofast - O3| / O3")
    ax[1].set_title(f"pixels changed by -Ofast ({len(xs)} of {a.size})")
    for x in ax:
        x.set_xlabel("α [M]"), x.set_ylabel("β [M]"), x.grid(False)
    save(fig, "raptor.png")


# ---------------------------------------------------------------- Q1
def parse_micro(txt, which):
    block = txt.split(which, 1)[1]
    rows = {}
    for line in block.splitlines()[3:9]:
        m = re.match(r"(\S+(?: double)?)\s+\d+ B \|(.*)\|(.*)\|", line)
        if m:
            rows[m.group(1)] = ([float(x) for x in m.group(2).split()], [float(x) for x in m.group(3).split()])
    return rows


def fig_q1():
    txt = read("q1_microbench.txt")
    base = parse_micro(txt, "Build: g++ -O2")
    nat = parse_micro(txt, "Build: g++ -O3 -march=native")
    ops = ["add", "mul", "div", "sqrt", "sin"]
    types = ["float16_t", "float", "double", "long double", "__float128"]
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    panels = [("latency, dependent chain (-O2)", base, 0), ("throughput, arrays (-O2, SSE2)", base, 1),
              ("throughput, arrays (-O3 -march=native)", nat, 1)]
    w = 0.16
    for ax, (title, data, k) in zip(axs, panels):
        for j, t in enumerate(types):
            if t not in data:
                continue
            ax.bar(np.arange(5) + (j - 2) * w, data[t][k], w * 0.88, color=SERIES[j], label=t)
        ax.set_yscale("log")
        ax.set_xticks(range(5), ops)
        ax.set_title(title)
    axs[0].set_ylabel("ns per operation (log)")
    axs[0].legend(ncol=1, fontsize=7, loc="upper left")
    save(fig, "q1_micro.png")


def fig_validate():
    rows = []
    for line in read("validate_shadow.txt").splitlines():
        m = re.match(r"(float|double|long double|__float128)\s+(\w+)\s+(\S+)\s+\d+ \| rel\. edge error: mean (\S+)\s+max (\S+)", line)
        if m:
            rows.append((f"{m.group(1)} {m.group(2)} {m.group(3)}", float(m.group(5)), m.group(1)))
    color = {"float": SERIES[1], "double": SERIES[0], "long double": SERIES[2], "__float128": SERIES[3]}
    fig, ax = plt.subplots(figsize=(7, 3.2))
    y = np.arange(len(rows))
    ax.barh(y, [r[1] for r in rows], color=[color[r[2]] for r in rows], height=0.7)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("max relative error of the shadow edge vs analytic curve (log)")
    for yi, r in zip(y, rows):
        ax.text(r[1] * 1.3, yi, f"{r[1]:.1e}", va="center", fontsize=7, color=INK2)
    save(fig, "validate.png")


# ---------------------------------------------------------------- parallel
def fig_parallel():
    txt = read("parallel.txt")
    blocks = re.split(r"\n(?=double, |float, )", txt)
    fig, ax = plt.subplots(figsize=(8, 4))
    data = {}
    for blk in blocks:
        m = re.match(r"(double|float), ", blk)
        if not m:
            continue
        rows = re.findall(r"^(\S.*?)\s{2,}([\d.]+)\s+([\d.]+)\s+(\d+)%", blk, re.M)
        data[m.group(1)] = [(r[0].strip(), float(r[2])) for r in rows]
    names = [n for n, _ in data["double"]]
    y = np.arange(len(names))
    h = 0.38
    for j, t in enumerate(["double", "float"]):
        if t in data:
            vals = dict(data[t])
            ax.barh(y + (j - 0.5) * h, [vals.get(n, 0) for n in names], h * 0.9, color=SERIES[j], label=t)
    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.axvline(4, color=INK2, lw=1, ls="--")
    ax.text(4.1, len(names) - 0.6, "4 cores", color=INK2, fontsize=7)
    ax.set_xlabel("speed-up over serial scalar double/float (same image)")
    ax.legend(loc="lower right")
    save(fig, "parallel.png")


def fig_cache():
    rows = re.findall(r"^\s*(\d+) (KiB|MiB) \|\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+) \|\s+([\d.]+)", read("cache_sweep.txt"), re.M)
    kib = np.array([int(r[0]) * (1024 if r[1] == "MiB" else 1) for r in rows])
    gf, gd = np.array([float(r[4]) for r in rows]), np.array([float(r[5]) for r in rows])
    lat = np.array([float(r[6]) for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.5))
    ax[0].plot(kib, gf, color=SERIES[1], marker="o", ms=4, label="float")
    ax[0].plot(kib, gd, color=SERIES[0], marker="o", ms=4, label="double")
    ax[0].set_ylabel("elements summed per ns")
    ax[0].set_title("streaming sum: float moves half the bytes")
    ax[1].plot(kib, lat, color=SERIES[2], marker="o", ms=4)
    ax[1].set_yscale("log")
    ax[1].set_ylabel("ns per dependent load (log)")
    ax[1].set_title("pointer chasing: load-to-use latency")
    for a in ax:
        a.set_xscale("log", base=2)
        a.set_xlabel("working set [KiB] (log2)")
        for size, lab in [(48, "L1d 48K"), (2048, "L2 2M"), (260 * 1024, "L3 260M")]:
            a.axvline(size, color=INK2, lw=0.8, ls=":")
            a.text(size * 1.08, a.get_ylim()[1] * 0.9 if a.get_yscale() == "linear" else a.get_ylim()[1] * 0.6, lab,
                   fontsize=7, color=INK2)
    ax[0].legend(loc="lower left")
    save(fig, "cache.png")


def fig_chaos():
    rows = re.findall(r"^\s+(\de-\d+)\s+(\d+)\s+(\S+)\s+(\S+)", read("q6_extra.txt"), re.M)
    d = np.array([float(r[0]) for r in rows])
    amp = np.array([float(r[3]) for r in rows])
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    ax.loglog(d, amp, color=SERIES[0], marker="o", ms=5, label="measured (double, DP45 rtol 1e-12)")
    ax.loglog(d, amp[0] * d[0] / d, color=INK2, lw=1, ls="--", label="∝ 1/δ")
    ax.invert_xaxis()
    ax.set_xlabel("δ = distance to the shadow edge (relative, log)")
    ax.set_ylabel("error amplification (log)")
    ax.legend(fontsize=7)
    save(fig, "chaos.png")


def fig_fastmath():
    a = kimg.load(os.path.join(BUILD, "fm_O3_float.bin"))
    b = kimg.load(os.path.join(BUILD, "fm_fast_float.bin"))
    rel = np.where(a["I"] > 0, np.abs(a["I"] - b["I"]) / np.maximum(a["I"], 1e-300), np.nan)
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    im = ax.imshow(np.log10(np.where(rel > 0, rel, np.nan)), origin="lower", cmap="Oranges",
                   extent=[-20, 20, -20, 20], vmin=-9, vmax=-4)
    fig.colorbar(im, ax=ax, fraction=0.046, label="log10 relative change")
    ax.set_title("float render: -O3 vs -O3 -ffast-math")
    ax.set_xlabel("α [M]"), ax.set_ylabel("β [M]"), ax.grid(False)
    save(fig, "fastmath.png")


if __name__ == "__main__":
    for f in [fig_render, fig_raptor, fig_q1, fig_validate, fig_parallel, fig_cache, fig_chaos, fig_fastmath]:
        try:
            f()
        except Exception as e:  # keep going so partial results still produce figures
            print(f"{f.__name__} failed: {e!r}")
