import sys, json, csv; sys.path.insert(0, ".")
from digit import *
SRC = "../fig/p27-1.png"
g, ink = load(SRC)
ax, bx, rx = calib([172.5, 401.5, 630.5, 861.5, 1088.0], [0, 200, 400, 600, 800])
ay, by, ry = calib([43.0, 157.5, 273.0, 388.5, 503.5, 619.5, 734.5, 850.5, 964.5], [16, 14, 12, 10, 8, 6, 4, 2, 0])
px = lambda re, xr: ((re - bx) / ax, (xr - by) / ay)
X0, X1, Y0, Y1 = 172, 1087, 41, 963
inside = np.zeros_like(ink); inside[Y0 + 4:Y1 - 3, X0 + 4:X1 - 3] = True
legend = np.zeros_like(ink); legend[:445, :910] = True
ticks_band = np.zeros_like(ink); ticks_band[Y0:Y1, X0:X0 + 18] = True; ticks_band[Y1 - 18:Y1, X0:X1] = True; ticks_band[Y0:Y1, X1 - 18:X1] = True
leg0 = np.zeros_like(ink); leg0[:445, :580] = True
pk = [refine_circle(ink, x, y, 11) for x, y, c in markers(g, (290, 422), 11, thr=0.6, mask=np.pad(np.ones((Y1 - Y0, X1 - X0 + 10), bool), ((Y0, ink.shape[0] - Y1), (X0, ink.shape[1] - X1 - 10))) & ~leg0)]
work = ink & inside & ~leg0
for cx, cy in pk: work[int(cy) - 13:int(cy) + 14, int(cx) - 13:int(cx) + 14] = False
SER = {
    "central": [(100, 3.15), (150, 4.13), (200, 5.0), (250, 5.84), (300, 6.68), (350, 7.39), (400, 7.95), (500, 9.25), (600, 10.5), (800, 12.8)],
    "hybrid":  [(100, 3.2), (300, 6.7), (400, 8.05), (600, 8.45), (800, 7.65)],
    "INS3D":   [(390, 8.1), (500, 9.75), (600, 11.0), (800, 13.45)],
    "KimMoin": [(100, 3.35), (150, 4.31), (200, 5.27), (250, 6.22), (300, 7.12), (350, 7.89), "px", (567, 512), (643, 457), (700, 430), (760, 397), (877, 343), "data", (600, 10.5), (800, 12.1)],
    "TEACH":   [(400, 8.1), (500, 8.55), (600, 8.2), (800, 7.45)],
}
out = {}; paths = []
for name, way in SER.items():
    gp, mode = [], "data"
    for w in way:
        if isinstance(w, str): mode = w; continue
        gp.append(px(*w) if mode == "data" else w)
    xs, ys = trace_fn(work, gp, gp[0][0], gp[-1][0], tol=6, blend=0.8)
    path = np.stack([xs, ys], 1); found = np.ones(len(xs), bool)
    span = gp[-1][0] - gp[0][0]
    data = [(ax * x + bx, ay * y + by) for x, y in path]
    out[name] = dict(points=len(xs), Re=[round(data[0][0]), round(data[-1][0])], columns_with_ink=f"{len(xs) / span:.0%}")
    paths.append(path)
    with open(f"fig10_{name}.csv", "w") as f:
        w = csv.writer(f); w.writerow(["Re", "Xr_over_h", "on_ink"]); [w.writerow([round(a, 2), round(b, 4), int(c)]) for (a, b), c in zip(data, found)]
with open("fig10_experiment.csv", "w") as f:
    w = csv.writer(f); w.writerow(["Re", "Xr_over_h"]); [w.writerow([round(ax * x + bx, 1), round(ay * y + by, 3)]) for x, y in pk]
cov, miss = coverage(work, paths, inside & ~legend & ~ticks_band)
print(json.dumps(out, indent=0)); print(f"curve ink explained: {cov:.1%}; calibration rms x {rx:.2f}px y {ry:.2f}px")
im = Image.open(SRC).convert("RGB"); d = ImageDraw.Draw(im)
cols = dict(central=(220, 0, 0), hybrid=(0, 110, 230), INS3D=(0, 160, 0), KimMoin=(230, 140, 0), TEACH=(150, 0, 200))
for (name, _), p in zip(SER.items(), paths): [d.point((x, y), fill=cols[name]) for x, y in p] or [d.ellipse([x - 1, y - 1, x + 1, y + 1], fill=cols[name]) for x, y in p]
for cx, cy in pk: d.ellipse([cx - 6, cy - 6, cx + 6, cy + 6], outline=(0, 160, 160), width=3)
ys, xs = np.nonzero(miss); [d.point((x, y), fill=(255, 0, 255)) for x, y in zip(xs, ys)]
im.save("fig10_overlay.png")

for name, y in (("central", 84), ("hybrid", 182), ("INS3D", 274), ("KimMoin", 322), ("TEACH", 372)):
    row = ink[y - 4:y + 5, 205:362].any(0); r = runs(row, 1)
    on = [e - s for s, e in r]; off = [r[i + 1][0] - r[i][1] for i in range(len(r) - 1)]
    print(f"{name:8s} legend sample: on {on}  off {off}")
print("markers:", len(pk))
