import sys, json, csv; sys.path.insert(0, ".")
from digit import *
SRC = "fig8_deskew.png"
g, ink = load(SRC)
X0, X1, Y0, Y1 = 431, 1627, 39, 1142
ax, bx, rx = calib([431.25, 729.0, 1025.75, 1326.25, 1626.75], [0, 100, 200, 300, 400])
ay, by, ry = calib([39.0, 312.5, 592.75, 868.5, 1142.25], [0, -1, -2, -3, -4])
px = lambda it, lv: ((it - bx) / ax, (lv - by) / ay)
inside = np.zeros_like(ink); inside[Y0 + 3:Y1 - 2, X0 + 3:X1 - 2] = True
legend = np.zeros_like(ink); legend[60:230, 1270:1560] = True
ticks_band = np.zeros_like(ink); ticks_band[Y0:Y1, X0:X0 + 3] = True; ticks_band[Y0:Y1, X1 - 18:X1] = True; 
for t in (729.0, 1025.75, 1326.25): ticks_band[Y0:Y0 + 18, int(t) - 4:int(t) + 5] = True; ticks_band[Y1 - 18:Y1, int(t) - 4:int(t) + 5] = True
ticks_band[Y1 - 3:Y1 + 2, X0:X1] = True
work = ink & inside & ~legend & ~ticks_band
paths = {}; out = {}
U0, U1 = np.array(px(0.3, -0.87)), np.array(px(219, -4.0)); P0, P1 = np.array(px(2, -1.35)), np.array(px(178, -4.0))
paths["u"] = chain(work, X0 + 1, X1 - 20, Y0 + 3, Y1 - 1, U0, U1)
um = np.zeros_like(ink); um[paths["u"][:, 1].astype(int), paths["u"][:, 0].astype(int)] = True
cst = np.where(ndi.distance_transform_edt(~um) <= 5, 4.0, 1.0)
VIA = [P0, np.array(px(7, -1.83)), np.array(px(65, -2.62)), P1]
paths["p"] = np.vstack([chain(work, X0 + 1, X1 - 20, Y0 + 3, Y1 - 1, a, b, cost=cst, gapw=8.0) for a, b in zip(VIA[:-1], VIA[1:])])
u = paths["u"]; pp = paths["p"]
taken = np.zeros_like(ink)
for p in (u, pp): taken[np.clip(p[:, 1].round().astype(int), 0, ink.shape[0] - 1), np.clip(p[:, 0].round().astype(int), 0, ink.shape[1] - 1)] = True
near = ndi.distance_transform_edt(~taken)
rest = work & (near > 3.5)
vy = np.full(len(u), np.nan)
for i, (x, y) in enumerate(u):
    col = np.nonzero(rest[int(max(y - 20, 0)):int(y + 21), int(x)])[0]
    if len(col): vy[i] = y - 20 + np.median(col)
ok = np.isfinite(vy); lab_, n_ = ndi.label(ok)
for k in range(1, n_ + 1):
    if (lab_ == k).sum() < 4: ok[lab_ == k] = False
v = u[:, 1].copy(); ii = np.nonzero(ok)[0]
off = np.zeros(len(u)); off[ok] = vy[ok] - u[ok, 1]
for a, b in zip(ii[:-1], ii[1:]):
    if 1 < b - a <= 30: off[a:b] = np.interp(np.arange(a, b), [a, b], [off[a], off[b]])
off = ndi.uniform_filter1d(off, 5); paths["v"] = np.stack([u[:, 0], u[:, 1] + off], 1)
for name in ("u", "v", "p"):
    p = paths[name]; data = [(ax * x + bx, ay * y + by) for x, y in p]
    with open(f"fig8_{name}.csv", "w") as f:
        w = csv.writer(f); w.writerow(["iteration", "log10_correction"]); [w.writerow([round(a, 2), round(b, 4)]) for a, b in data]
    out[name] = dict(points=len(p), iterations=[round(data[0][0], 1), round(data[-1][0], 1)], start=round(data[0][1], 3))
core = inside & ~legend & ~ticks_band
cov, miss = coverage(work, list(paths.values()), core)
print(json.dumps(out, indent=0)); print(f"curve ink explained {cov:.1%}; calibration rms x {rx:.2f}px y {ry:.2f}px")
im = Image.open(SRC).convert("RGB"); d = ImageDraw.Draw(im)
for (name, p), c in zip(paths.items(), ((220, 0, 0), (0, 160, 0), (0, 110, 230))): d.line([tuple(q) for q in p], fill=c, width=2)
ys, xs = np.nonzero(miss); [d.point((x, y), fill=(255, 0, 255)) for x, y in zip(xs, ys)]
im.save("fig8_overlay.png")
