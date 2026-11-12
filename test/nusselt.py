import sys, numpy as np
L, M = 161, 21; re = float(sys.argv[2])
a = np.array([[float(l[0:8]), float(l[8:16]), float(l[24:36]), float(l[73:84])] for l in open(sys.argv[1]) if len(l) > 100]).reshape(4, L, M, 4)[1]
x = a[:, 0, 0] if len(sys.argv) < 4 else a[:, 0, 0] + 0.5 * float(sys.argv[3]); y = a[0, :, 1]; u = a[:, :, 2]; t = a[:, :, 3]
tb = np.trapezoid(u * t, y, axis=1) / np.trapezoid(u, y, axis=1)
g = (-3 * t[:, 0] + 4 * t[:, 1] - t[:, 2]) / (2 * (y[1] - y[0]))
nu_w = 2 * g / tb
sel = (x > 2) & (x < 6)
lam = -np.polyfit(x[sel], np.log(tb[sel]), 1)[0]
print(f"Re {re:g}: centre/mean u {u[L // 2, M // 2] / np.trapezoid(u[L // 2], y):.3f}, Nu from the bulk decay {lam * re:.2f}, "
      f"from the wall gradient {np.mean(nu_w[sel]):.2f} (exact 7.54), Tb at x = 2, 4, 6: {np.interp([2, 4, 6], x, tb).round(4).tolist()}")
