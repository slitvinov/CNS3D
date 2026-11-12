import numpy as np, sys
a = sys.argv[1:]; M, b = (int(a[0]), float(a[1])) if len(a) > 1 else (43, 2.5)
qw, ql, bu, NU, LD, tu, ls = (float(a[2]), float(a[3]), float(a[4]), int(a[5]), float(a[6]), float(a[7]), float(a[8])) if len(a) > 8 else (None, None, 0.0, 9, 30.0, 0.04, 0.14)
L, NZ = 45, 4
su = np.linspace(0, 1, NU); xu = -4 * (np.exp(bu * (1 - su)) - 1) / (np.exp(bu) - 1) if bu > 0 else np.linspace(-4, 0, NU)
s = np.linspace(0, 1, L - NU + 1); xd = LD * (np.exp(b * s) - 1) / (np.exp(b) - 1)
x = np.concatenate([xu, xd[1:]]); dz = 3.0 / (M - 1); Is = NU
def seg(y0, y1, n, qa, qb):
    t = (np.arange(n) + 0.5) / n; w = 1 / (1 + qa * np.exp(-8 * t) + qb * np.exp(-8 * (1 - t)))
    return y0 + (y1 - y0) * np.concatenate([[0], np.cumsum(w) / w.sum()])
Js = 15
y = np.concatenate([seg(0.0, 1.0, Js - 1, qw or 0.0, ql or 0.0), seg(1.0, 3.0, M - Js, ql or 0.0, qw or 0.0)[1:]])
k0 = 1.5 * tu ** 2; e0 = 0.09 ** 0.75 * k0 ** 1.5 / ls
o = [f"{L} {M} {NZ} 1 {Is} 1 {Js} 1 {NZ}"]
for k in range(NZ):
    for i in range(L):
        for j in range(M):
            wall = j == 0 or j == M - 1 or (i < Is and j < Js)
            u = 0.0 if wall else (1.0 if x[i] <= 0 else 2 / 3)
            o.append(f"{x[i]:.6f} {y[j]:.6f} {(k - 0.5) * dz:.6f} {u:.6e} 0 0 {k0:.4e} {e0:.4e}")
open("field.txt", "w").write("\n".join(o) + "\n"); print("step corner I,J =", Is, Js, "y=", y[Js - 1])
