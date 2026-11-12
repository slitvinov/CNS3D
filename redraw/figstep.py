import numpy as np, json, os
def load(f, L, M):
    r = [l for l in open(f) if len(l) > 100]
    a = np.array([[float(l[0:8]), float(l[8:16]), float(l[24:36]), float(l[36:48]), float(l[84:95])] for l in r]).reshape(4, L, M, 5)[1]
    x = a[:, 0, 0]; y = a[0, :, 1]; u = a[:, :, 2]; k = a[:, :, 4]
    psi = np.zeros((L, M)); psi[:, 1:] = np.cumsum(0.5 * (u[:, 1:] + u[:, :-1]) * np.diff(y)[None, :], axis=1)
    return x, y, u, psi / psi[-1, -1], k
def cut(x, y, F, x0, x1):
    i = (x > x0) & (x < x1); xs = np.r_[x0, x[i], x1]
    edge = lambda xe: [np.interp(xe, x, F[:, j]) for j in range(len(y))]
    G = np.vstack([edge(x0), F[i], edge(x1)])
    return xs, G
def write(fn, xs, y, G):
    with open(fn, "w") as f:
        f.write(f"{len(xs)} {len(y)}\n" + " ".join(f"{v:.5f}" for v in xs) + "\n" + " ".join(f"{v:.5f}" for v in y) + "\n")
        for j in range(len(y)): f.write(" ".join(f"{G[i, j]:.6e}" for i in range(len(xs))) + "\n")
def xr(x, u):
    i = next((i for i in range(2, len(x)) if u[i - 1, 1] < 0 <= u[i, 1]), None)
    return float(np.interp(0, [u[i - 1, 1], u[i, 1]], [x[i - 1], x[i]]))
out = {}
for tag, f in (("a", "step_hyb600.f7"), ("b", "step_pc20600.f7")):
    if not os.path.exists(f): continue
    x, y, u, psi, k = load(f, 45, 45); xs, G = cut(x, y, psi, 0.0, 23.0); write(f"fig11{tag}.dat", xs, y, G)
    out[f"fig11{tag}"] = dict(Xr=round(xr(x, u), 2), psi_min=round(float(psi.min()), 4), psi_max=round(float(psi.max()), 4))
x, y, u, psi, k = load("tstep.f7", 45, 42)
for nm, F_ in (("fig14", psi), ("fig15", k)):
    xs, G = cut(x, y, F_, -4.0, 0.0); write(nm + "u.dat", xs, y, G)
    xs, G = cut(x, y, F_, 0.0, 20.0); write(nm + "d.dat", xs, y, G)
i0 = int(np.argmin(abs(x))); out["fig14"] = dict(Xr=round(xr(x[i0:], u[i0:]), 2), psi_min=round(float(psi.min()), 4)); out["fig15"] = dict(k_max=round(float(k[x > 0].max()), 4))
json.dump(out, open("figstep.json", "w")); print(out)
