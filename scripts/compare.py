#!/usr/bin/env python3
"""Small diff utilities used by run_all.sh.

  compare.py libm  A.bin B.bin [labelA labelB]   per-function mismatch counts (q5_env output)
  compare.py kimg  A.bin B.bin [labelA labelB]   pixel-level diff of two KIMG renders
  compare.py raptor A.dat B.dat                  diff of two RAPTOR img_data files
"""
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import kimg  # noqa: E402

FUNCS = ["sin", "cos", "tan", "exp", "log", "pow", "atan2", "cbrt", "sinh", "lgamma"]
N = 1_000_000


def ulp_diff(a, b):
    """distance in units of the last place between two float64 arrays"""
    ia = a.view(np.int64).astype(np.int64)
    ib = b.view(np.int64).astype(np.int64)
    ia = np.where(ia < 0, np.int64(-(2**63)) - ia, ia)
    ib = np.where(ib < 0, np.int64(-(2**63)) - ib, ib)
    return np.abs(ia - ib)


def libm(fa, fb, la="A", lb="B"):
    a = np.fromfile(fa, dtype=np.float64).reshape(len(FUNCS), N)
    b = np.fromfile(fb, dtype=np.float64).reshape(len(FUNCS), N)
    print(f"libm results that differ between {la} and {lb} ({N} inputs each):")
    for k, name in enumerate(FUNCS):
        both = np.isfinite(a[k]) & np.isfinite(b[k])
        d = a[k][both] != b[k][both]
        u = ulp_diff(a[k][both], b[k][both])
        print(f"  {name:7s} {int(d.sum()):8d} differ ({100*d.mean():6.3f}%)  max {int(u.max()) if u.size else 0} ulp")


def kimgdiff(fa, fb, la="A", lb="B"):
    a, b = kimg.load(fa), kimg.load(fb)
    n = a["I"].size
    st = int((a["status"] != b["status"]).sum())
    bit = int((a["I"] != b["I"]).sum())
    both = (a["status"] == 3) & (b["status"] == 3)
    rel = np.abs(a["I"] - b["I"])[both] / np.maximum(np.abs(a["I"][both]), 1e-300)
    th = int((a["theta"] != b["theta"]).sum())
    print(f"{la} vs {lb}: {n} pixels, status differs {st}, intensity bitwise differs {bit} "
          f"({100*bit/n:.1f}%), final theta differs {th}, max rel. intensity diff {rel.max():.3e}, "
          f"flux {a['I'].sum():.17g} vs {b['I'].sum():.17g}")


def raptor(fa, fb):
    def load(f):
        with open(f) as h:
            head = h.readline().split()
            v = np.array([float(x) for x in h.read().split()])
        return head, v
    ha, a = load(fa)
    hb, b = load(fb)
    d = a != b
    rel = np.abs(a - b)[a != 0] / np.abs(a[a != 0])
    print(f"RAPTOR image {ha[0]}x{ha[1]}: total flux {ha[2]} Jy vs {hb[2]} Jy")
    print(f"  pixels differing in the printed 6 significant digits: {int(d.sum())} of {a.size} "
          f"({100*d.mean():.2f}%), max rel. diff {rel.max():.3e}, median rel. diff of differing pixels "
          f"{np.median(rel[rel > 0]) if (rel > 0).any() else 0:.3e}")


if __name__ == "__main__":
    cmd, *args = sys.argv[1:]
    {"libm": libm, "kimg": kimgdiff, "raptor": raptor}[cmd](*args)
