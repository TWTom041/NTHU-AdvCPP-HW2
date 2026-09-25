"""Reader for the KIMG binary images written by kerr::Image::save."""
import numpy as np


def load(path):
    with open(path, "rb") as f:
        assert f.read(4) == b"KIMG"
        w, h = np.frombuffer(f.read(8), dtype=np.int32)
        n = int(w) * int(h)
        inten = np.frombuffer(f.read(8 * n), dtype=np.float64).reshape(h, w)
        status = np.frombuffer(f.read(n), dtype=np.uint8).reshape(h, w)
        theta = np.frombuffer(f.read(8 * n), dtype=np.float64).reshape(h, w)
        phi = np.frombuffer(f.read(8 * n), dtype=np.float64).reshape(h, w)
        steps = np.frombuffer(f.read(8 * n), dtype=np.int64).reshape(h, w)
    return dict(w=int(w), h=int(h), I=inten, status=status, theta=theta, phi=phi, steps=steps)


STATUS = ["running", "captured", "escaped", "disk", "maxsteps", "notfinite"]
