import numpy as np, sys, os
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "redraw")
COL = ((0, 8), (8, 8), (16, 8), (24, 12), (36, 12), (48, 12), (62, 11), (84, 11), (95, 11))
STATIONS = (0, 60, 90)
RC = 2.3
S_UP, S_DN = -1.031, RC * np.pi / 2 + 1.0906
PLANES = (("a", 0.0), ("b", 0.2451))

def read(d, suf):
    h = open(f"{d}/field").readline().split()
    L, M, N = map(int, h[:3])
    a = np.array([[float(l[p:p + w]) for p, w in COL] for l in open(f"{d}/field{suf}") if len(l) > 100]).reshape(N, L, M, 9)
    g = np.array([[float(l[p:p + w]) for p, w in COL[:3]] for l in open(f"{d}/field") if len(l) > 100]).reshape(N, L, M, 3)
    a[..., :3] = g
    return a, L, M, N, h

def arc(a):
    M = a.shape[2]
    x, y = a[0, :, M // 2, 0], a[0, :, M // 2, 1]
    return np.where(x < 0, x, np.where(y > 0, RC * np.pi / 2 + y, RC * np.arctan2(x, -y)))

def at(a, s0):
    s = arc(a)
    i = int(np.abs(s - s0).argmin())
    if abs(s[i] - s0) > 2e-3:
        sys.exit(f"no section at {s0:.4f}")
    return i

def amean_at(a, s0, c):
    s = arc(a)
    i = min(max(int(np.searchsorted(s, s0)), 1), len(s) - 1)
    f = (s0 - s[i - 1]) / (s[i] - s[i - 1])
    return amean(a, i - 1, c) * (1 - f) + amean(a, i, c) * f

def plane(a, i, zz):
    z = a[:, i, 1, 2]
    k = min(max(int(np.searchsorted(z, zz)), 1), len(z) - 1)
    f = (zz - z[k - 1]) / (z[k] - z[k - 1])
    return a[k - 1, i] * (1 - f) + a[k, i] * f

def amean(a, i, c):
    s = a[:, i]
    N, M = s.shape[:2]
    sj = np.concatenate([[0], np.cumsum(np.hypot(np.diff(s[0, :, 0]), np.diff(s[0, :, 1])))])
    z = s[:, 0, 2]
    wj = np.zeros(M)
    wj[1:M - 1] = (sj[2:] - sj[:-2]) / 2
    wj[1] += (sj[1] - sj[0]) / 2
    wj[M - 2] += (sj[M - 1] - sj[M - 2]) / 2
    wk = np.zeros(N)
    wk[0] = (z[1] - z[0]) / 2
    wk[1:N - 1] = (z[2:] - z[:-2]) / 2
    wk[N - 2] += (z[N - 1] - z[N - 2]) / 2
    W = np.outer(wk, wj)
    return (W * s[..., c]).sum() / W.sum()

def profiles(a):
    out = {}
    for th in STATIONS:
        t = np.radians(th)
        i = at(a, RC * t)
        for pn, zz in PLANES:
            s = plane(a, i, zz)
            us = s[:, 3] * np.cos(t) + s[:, 4] * np.sin(t)
            rd = np.hypot(s[:, 0], s[:, 1]) - 1.8
            o = np.argsort(rd)
            out[(th, pn)] = (rd[o], us[o])
    return out

def bend(d, suf):
    a, L, M, N, h = read(d, suf)
    P = profiles(a)
    rms = {}
    turb = open(f"{d}/control").readline()[5:10].strip() == "2"
    for src in () if turb else ("present", "humphrey", "rhie"):
        dev = []
        for (th, pn), (rd, us) in P.items():
            q = np.loadtxt(f"{R}/fig19{pn}_theta{th}_{src}.csv", delimiter=",", skiprows=1)
            dev.append(np.interp(q[:, 1], rd, us) - q[:, 0])
        rms[src] = np.sqrt(np.mean(np.concatenate(dev) ** 2))
    dp = amean_at(a, S_UP, 6) - amean_at(a, S_DN, 6)
    db = amean(a, at(a, 0), 6) - amean(a, at(a, RC * np.pi / 2), 6)
    s = a[:, at(a, RC * np.pi / 2)]
    ur = np.abs(s[..., 3]).max()
    w = np.abs(s[..., 5]).max()
    k = a[..., 7]
    yp = []
    for j, jw in ((1, 0), (M - 2, M - 1)):
        dist = np.hypot(a[:N - 1, :, j, 0] - a[:N - 1, :, jw, 0], a[:N - 1, :, j, 1] - a[:N - 1, :, jw, 1])
        yp.append((dist * np.sqrt(np.maximum(k[:N - 1, :, j], 0))).ravel())
    dist = np.abs(a[N - 2, :, 1:M - 1, 2] - a[N - 1, :, 1:M - 1, 2])
    yp.append((dist * np.sqrt(np.maximum(k[N - 2, :, 1:M - 1], 0))).ravel())
    re = float(open(f"{d}/control").readlines()[2][:7])
    yk = np.median(np.concatenate(yp)) * 0.09 ** 0.25 * re
    print(f"{d} dp {dp:.4f} db {db:.4f} ur {ur:.3f} w {w:.3f}" + (f" y+ {yk:.1f}" if turb else f" present {rms['present']:.3f} humphrey {rms['humphrey']:.3f} rhie {rms['rhie']:.3f}"))

def xr(d, suf=".o"):
    a, L, M, N, h = read(d, suf)
    x = a[1, :, 0, 0]
    u = a[1, :, 1, 3]
    z = [float(np.interp(0, [u[i - 1], u[i]], [x[i - 1], x[i]])) for i in range(max(int(h[4]), 2), L) if u[i - 1] < 0 <= u[i]]
    return z[0] if z else None


def step(d, suf):
    z = xr(d, suf)
    print(f"{d} xr {z:.3f}" if z is not None else f"{d} xr none")

def diff(ref, ds, suf):
    a, L, M, N, h = read(ref, suf)
    P0 = profiles(a)
    x = np.linspace(0.02, 0.98, 49)
    for d in ds:
        a, L, M, N, h = read(d, suf)
        P = profiles(a)
        print(d, " ".join(f"{th}{pn} {np.sqrt(np.mean((np.interp(x, *P[(th, pn)]) - np.interp(x, *P0[(th, pn)])) ** 2)):.3f}"
                          for th in STATIONS for pn, _ in PLANES))

if __name__ == "__main__":
    a = sys.argv[1:]
    suf = ".o"
    if a and a[0] == "-s":
        suf, a = a[1], a[2:]
    if a and a[0] == "-d":
        diff(a[1], a[2:], suf)
    else:
        for d in a:
            if int(open(f"{d}/field").readline().split()[2]) > 4:
                bend(d, suf)
            else:
                step(d, suf)
