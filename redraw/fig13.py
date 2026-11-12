import sys, json, csv; sys.path.insert(0, ".")
from digit import *
SRC = "fig13_deskew.png"
g, ink = load(SRC)
X0, Y0 = 157, 866
ax, bx, rx = calib([156.0, 423.0, 688.0, 945.0, 1207.0], [0, 2, 4, 6, 8])
ay, by, ry = calib([865.5, 734.5, 604.0, 471.0, 338.0, 205.5], [0, 0.2, 0.4, 0.6, 0.8, 1.0])
px = lambda x, y: ((x - bx) / ax, (y - by) / ay)
leg = np.zeros_like(ink); leg[:120, :] = True
inset = np.zeros_like(ink); inset[150:410, 500:1226] = True
sym = refine_circle(ink, 294, 80, 11)
pk = [refine_circle(ink, x, y, 11) for x, y, c in markers(g, (int(round(sym[0])), int(round(sym[1]))), 11, thr=0.5, mask=~leg & ~inset)
      if x > X0 + 10 and y < Y0 - 10]
work = ink.copy(); work[leg | inset] = False; work[:, :X0 + 6] = False; work[Y0 - 3:, :] = False
for cx, cy in pk: work[int(cy) - 13:int(cy) + 14, int(cx) - 13:int(cx) + 14] = False
lab_, n_ = ndi.label(work); big = np.argmax(ndi.sum(work, lab_, range(1, n_ + 1))) + 1
cm = lab_ == big; ys, xs = np.nonzero(cm)
a = np.array([xs.min(), ys[xs == xs.min()].mean()]); b = np.array([xs[ys == ys.max()].mean(), ys.max()])
path = chain(cm, xs.min() - 2, xs.max() + 3, ys.min() - 2, ys.max() + 3, a, b)
data = [(ax * x + bx, ay * y + by) for x, y in path]
with open("fig13_present.csv", "w") as f:
    w = csv.writer(f); w.writerow(["x_over_h", "y_over_h"]); [w.writerow([round(u, 4), round(v, 4)]) for u, v in data]
pts = sorted((ax * x + bx, ay * y + by) for x, y in pk)
with open("fig13_exp.csv", "w") as f:
    w = csv.writer(f); w.writerow(["x_over_h", "y_over_h"]); [w.writerow([round(u, 4), round(v, 4)]) for u, v in pts]
core = np.zeros_like(ink); core[120:Y0 - 3, X0 + 6:] = True; core &= ~inset
cov, miss = coverage(work, [path], core)
print(json.dumps(dict(curve=dict(points=len(path), start=[round(v, 3) for v in data[0]], end=[round(v, 3) for v in data[-1]]),
      exp=[[round(u, 3), round(v, 3)] for u, v in pts], legend_symbol=[round(v, 1) for v in sym])))
print(f"curve ink explained {cov:.1%}; calibration rms x {rx:.2f}px y {ry:.2f}px")
im = Image.open(SRC).convert("RGB"); d = ImageDraw.Draw(im)
d.line([tuple(q) for q in path], fill=(220, 0, 0), width=2)
for cx, cy in pk: d.ellipse([cx - 6, cy - 6, cx + 6, cy + 6], outline=(0, 160, 160), width=3)
ys, xs = np.nonzero(miss); [d.point((x, y), fill=(255, 0, 255)) for x, y in zip(xs, ys)]
im.save("fig13_overlay.png")
