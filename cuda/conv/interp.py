import numpy as np, sys

COL = ((24, 12), (36, 12), (48, 12), (62, 11), (73, 11), (84, 11), (95, 11))


def size(f):
    return tuple(map(int, open(f).readline().split()[:3]))


def values(f, n):
    v = np.array([[float(l[p:p + w]) for p, w in COL] for l in open(f) if len(l) > 100])
    return v.reshape(n[2], n[0], n[1], len(COL))


def coords(f, n):
    g = np.array([[float(l[p:p + 8]) for p in (0, 8, 16)] for l in open(f) if len(l) > 100]).reshape(n[2], n[0], n[1], 3)
    return g[0, :, 0, 0], g[0, 0, :, 1], g[:, 0, 0, 2]


def axis(nf, nc, xf=None, xc=None):
    if xf is not None:
        i = np.clip(np.searchsorted(xc, xf) - 1, 0, nc - 2)
        return i, np.clip((xf - xc[i]) / (xc[i + 1] - xc[i]), 0, 1)
    if (nf - 1) % (nc - 1):
        sys.exit(f"{nf} points are not a refinement of {nc}")
    t = np.arange(nf) * (nc - 1) / (nf - 1)
    i = np.minimum(t.astype(int), nc - 2)
    return i, t - i


phys = sys.argv[1] == "-x"
coarse, fine = sys.argv[1 + phys:3 + phys]
nc, nf = size(f"{coarse}/field"), size(f"{fine}/field")
a = values(f"{coarse}/field.o", nc)
cc, cf = (coords(f"{coarse}/field", nc), coords(f"{fine}/field", nf)) if phys else ((None,) * 3, (None,) * 3)
ii, fi = axis(nf[0], nc[0], cf[0], cc[0])
jj, fj = axis(nf[1], nc[1], cf[1], cc[1])
kk, fk = (np.arange(nf[2]), np.zeros(nf[2])) if nf[2] == nc[2] else axis(nf[2], nc[2])
b = np.zeros((nf[2], nf[0], nf[1], len(COL)))
for dk in (0, 1):
    wk = fk if dk else 1 - fk
    if not wk.any():
        continue
    for di in (0, 1):
        wi = fi if di else 1 - fi
        for dj in (0, 1):
            wj = fj if dj else 1 - fj
            w = wk[:, None, None] * wi[None, :, None] * wj[None, None, :]
            b += w[..., None] * a[np.ix_(np.minimum(kk + dk, nc[2] - 1), ii + di, jj + dj)]
lines = open(f"{fine}/field").read().split("\n")
out, n = [], 0
for l in lines:
    if len(l) > 100:
        k, r = divmod(n, nf[0] * nf[1])
        i, j = divmod(r, nf[1])
        if i == 0:
            out.append(l)
        else:
            u = b[k, i, j]
            out.append(l[:24] + "%12.4E%12.4E%12.4E  " % tuple(u[:3]) + "%11.4E%11.4E%11.4E%11.4E" % tuple(u[3:]) + l[106:])
        n += 1
    else:
        out.append(l)
open(f"{fine}/field", "w").write("\n".join(out))
