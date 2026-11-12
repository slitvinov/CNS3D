import sys, json, csv; sys.path.insert(0, ".")
from digit import *
SRC = "fig7_deskew.png"
g, ink = load(SRC)
X0, X1, Y0, Y1 = 213, 1434, 43, 1067
ax, bx, rx = calib([213.5, 389.5, 564.0, 740.5, 912.0, 1086.5, 1260.5, 1435.0], [-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
ay, by, ry = calib([1067.0, 859.5, 656.5, 450.5, 244.5, 42.5], [0, 0.2, 0.4, 0.6, 0.8, 1.0])
px = lambda u, y: ((u - bx) / ax, (y - by) / ay)
inside = np.zeros_like(ink); inside[Y0 + 4:Y1 - 3, X0 + 4:X1 - 3] = True
legend = np.zeros_like(ink); legend[515:810, 725:1400] = True
mk = np.zeros_like(ink); mk[Y0 - 5:Y1 + 5, X0 - 5:X1 + 5] = True; mk &= ~legend
lx, ly = 797, 724
pk = [refine_circle(ink, x, y, 11) for x, y, c in markers(g, (lx, ly), 11, thr=0.45, mask=mk)]
work = ink & inside & ~legend
for cx, cy in pk: work[int(cy) - 13:int(cy) + 14, int(cx) - 13:int(cx) + 14] = False
SER = {
    "orthogonal":    [(0, 0), (-0.06, 0.05), (-0.2, 0.19), (-0.245, 0.3), (-0.22, 0.4), (-0.13, 0.5), (-0.02, 0.59), (0.06, 0.655), (0.16, 0.75), (0.25, 0.86), (0.31, 0.9)],
    "nonorthogonal": [(0, 0), (-0.06, 0.05), (-0.19, 0.19), (-0.226, 0.3), (-0.2, 0.4), (-0.125, 0.5), (-0.02, 0.59), (0.05, 0.655), (0.14, 0.75), (0.235, 0.86), (0.31, 0.9)],
}
TOP = [(0.31, 0.9), (0.4, 0.925), (0.5, 0.945), (0.7, 0.97), (0.9, 0.985), (0.995, 0.99)]
out, paths = {}, []
gO = sorted(px(*w)[::-1] for w in SER["orthogonal"]); gN = sorted(px(*w)[::-1] for w in SER["nonorthogonal"])
GO = lambda r: np.interp(r, [q[0] for q in gO], [q[1] for q in gO]); GN = lambda r: np.interp(r, [q[0] for q in gN], [q[1] for q in gN])
ycross = (0.6 - by) / ay
rows = range(int(gO[0][0]), int(gO[-1][0]) + 1)
res = {"orthogonal": [], "nonorthogonal": []}
for r in rows:
    rr = [((s0 + e0 - 1) / 2, e0 - s0) for s0, e0 in runs(ink[r], 1) if 2 <= e0 - s0 <= 16]
    rr = [(x, w) for x, w in rr if all(np.hypot(x - cx, r - cy) > 12.5 for cx, cy in pk)]
    c = sorted(x for x, w in rr if min(abs(x - GO(r)), abs(x - GN(r))) < 22)
    if not c: continue
    if len(c) == 1:
        x = c[0]
        if abs(GO(r) - GN(r)) < 7: res["orthogonal"].append((x, r)); res["nonorthogonal"].append((x, r))
        elif abs(x - GO(r)) < abs(x - GN(r)): res["orthogonal"].append((x, r))
        else: res["nonorthogonal"].append((x, r))
        continue
    lo, hi = c[0], c[-1]
    o, n = (lo, hi) if r > ycross else (hi, lo)
    res["orthogonal"].append((o, r)); res["nonorthogonal"].append((n, r))
xs2, ys2 = [], []
for seg in ([(0.31, 0.9), (0.4, 0.925), (0.5, 0.945), (0.6, 0.958), (0.68, 0.967)], [(0.72, 0.978), (0.76, 0.983), (0.85, 0.990), (0.92, 0.993)]):
    lo = [px(*w) for w in seg]; a_, b_ = trace_fn(work, lo, lo[0][0], lo[-1][0], tol=6, blend=0.8); xs2 += list(a_); ys2 += list(b_)
xs2, ys2 = np.array(xs2), np.array(ys2)
for name in res:
    P = np.array(sorted(res[name], key=lambda q: -q[1]))
    allr = np.arange(int(P[:, 1].max()), int(P[:, 1].min()) - 1, -1)
    xi = np.interp(-allr, -P[:, 1], P[:, 0]); flag = np.isin(allr, P[:, 1]).astype(int)
    path = np.vstack([np.stack([xi, allr], 1), np.stack([xs2, ys2], 1)]); fl = np.r_[flag, np.ones(len(xs2), int)]
    data = [(ax * x + bx, ay * y + by) for x, y in path]
    out[name] = dict(points=len(path), measured=f"{fl.mean():.0%}")
    paths.append(path[fl == 1])
    with open(f"fig7_{name}.csv", "w") as f:
        w = csv.writer(f); w.writerow(["u_over_U", "y_over_w", "on_ink"]); [w.writerow([round(a, 4), round(b, 4), k]) for (a, b), k in zip(data, fl)]
pts = sorted((ay * y + by, ax * x + bx) for x, y in pk)
with open("fig7_burggraf.csv", "w") as f:
    w = csv.writer(f); w.writerow(["u_over_U", "y_over_w"]); [w.writerow([round(u, 4), round(y, 4)]) for y, u in pts]
core = np.zeros_like(ink); core[Y0 + 18:Y1 - 18, X0 + 18:X1 - 18] = True
cov, miss = coverage(work, paths, core & ~legend)
print(json.dumps(out, indent=0)); print(f"{len(pk)} points; curve ink explained {cov:.1%}; calibration rms x {rx:.2f}px y {ry:.2f}px")
im = Image.open(SRC).convert("RGB"); d = ImageDraw.Draw(im)
for p, c in zip(paths, ((220, 0, 0), (0, 110, 230))): [d.point((x, y), fill=c) for x, y in p]
for cx, cy in pk: d.ellipse([cx - 6, cy - 6, cx + 6, cy + 6], outline=(0, 160, 160), width=3)
ys, xs = np.nonzero(miss); [d.point((x, y), fill=(255, 0, 255)) for x, y in zip(xs, ys)]
im.save("fig7_overlay.png")
