import sys, json, csv; sys.path.insert(0, ".")
from digit import *
from skimage.measure import label, regionprops
SRC = "fig19_deskew.png"
g, ink = load(SRC)
XL, XR, SEP = 241.5, 1585.5, (682.0, 1130.5)
PAN = {"a": (45.5, 461.0), "b": (711.0, 1130.5)}
UPX = (570.0 - 241.5) / 1.5
ORG = (XL, SEP[0], SEP[1])
leg = np.zeros_like(ink); leg[470:700, 700:1772] = True
cand = markers(g, (811, 654), 11, thr=0.5, mask=~leg)
pk = [refine_circle(ink, x, y, 11) for x, y, c in cand if XL - 5 < x < XR + 5 and any(a - 5 < y < b + 5 for a, b in PAN.values())]
work = ink.copy(); work[leg] = False
for cx, cy in pk: yy, xx = np.ogrid[:ink.shape[0], :ink.shape[1]]; work[(yy - cy) ** 2 + (xx - cx) ** 2 <= 10.5 ** 2] = False
for a, b in PAN.values():
    work[int(a) - 2:int(a) + 3, :] = False; work[int(b) - 2:int(b) + 3, :] = False
    work[int(a):int(b), int(XL) - 4:int(XL) + 16] = False; work[int(a):int(b), int(XR) - 16:int(XR) + 5] = False
    for x in SEP: work[int(a):int(a) + 16, int(x) - 3:int(x) + 4] = False; work[int(b) - 16:int(b), int(x) - 3:int(x) + 4] = False
work[:, :int(XL) - 4] = False; work[:, int(XR) + 5:] = False
for x0_, y0_, x1_, y1_ in ((385, 228, 505, 282), (755, 345, 895, 398), (1145, 332, 1285, 385),
                           (385, 893, 505, 948), (730, 893, 875, 952), (1115, 885, 1260, 942)):
    work[y0_:y1_, x0_:x1_] = False
curve = work
ringzone = np.zeros_like(ink); yy, xx = np.ogrid[:ink.shape[0], :ink.shape[1]]
for cx, cy in pk: ringzone |= (yy - cy) ** 2 + (xx - cx) ** 2 <= 14 ** 2
out = {}; paths = []
for pn, (ya, yb) in PAN.items():
    for j, (th, x0) in enumerate(zip((0, 60, 90), ORG)):
        x1 = ORG[j + 1] if j < 2 else XR
        pres = chain(curve, x0 - 2, x1 - 6, ya - 1, yb + 2, np.array([x0, ya]), np.array([x0, yb]))
        if len(pres) < 10: continue
        pm = np.zeros_like(ink); pm[np.clip(pres[:, 1].round().astype(int), 0, ink.shape[0] - 1), np.clip(pres[:, 0].round().astype(int), 0, ink.shape[1] - 1)] = True
        far = curve & (ndi.distance_transform_edt(~pm) > 4)
        far = ndi.binary_opening(far, iterations=1)
        far &= ~ringzone
        far[:, :int(x0) - 2] = False; far[:, int(x1) - 6:] = False; far[:int(ya), :] = False; far[int(yb) + 1:, :] = False
        fy, fx = np.nonzero(far)
        near_far = ndi.distance_transform_edt(~far)
        keep = near_far[np.clip(pres[:, 1].round().astype(int), 0, ink.shape[0] - 1), np.clip(pres[:, 0].round().astype(int), 0, ink.shape[1] - 1)] > 30
        sm = np.stack([ndi.uniform_filter1d(pres[:, k], 9, mode="nearest") for k in (0, 1)], 1)
        tg = np.gradient(sm, axis=0); tg /= np.linalg.norm(tg, axis=1)[:, None] + 1e-9; nm = np.stack([-tg[:, 1], tg[:, 0]], 1)
        fy, fx = np.nonzero(far); F = np.stack([fx, fy], 1).astype(float); off = np.full(len(pres), np.nan)
        if len(F):
            from scipy.spatial import cKDTree
            T = cKDTree(F)
            for i, (q, t, n) in enumerate(zip(sm, tg, nm)):
                idx = T.query_ball_point(q, 25)
                if not idx: continue
                v = F[idx] - q; along = v @ t; acr = v @ n; sel = np.abs(along) <= 3
                if sel.any(): off[i] = np.median(acr[sel])
        ok = np.isfinite(off)
        lab_, n_ = ndi.label(ok)
        for k_ in range(1, n_ + 1):
            if (lab_ == k_).sum() < 5: ok[lab_ == k_] = False
        off[~ok] = np.nan
        if ok.any(): off[ok] = ndi.median_filter(np.where(ok, off, 0), 7)[ok]
        o = np.zeros(len(pres))
        if ok.any():
            ii = np.nonzero(ok)[0]; o[ok] = off[ok]
            for a_, b_ in zip(ii[:-1], ii[1:]):
                if 1 < b_ - a_ <= 40: o[a_:b_] = np.interp(np.arange(a_, b_), [a_, b_], [off[a_], off[b_]])
        o = ndi.uniform_filter1d(o, 5)
        rhie = sm + nm * o[:, None]
        for kind, p in (("present", pres), ("rhie", rhie)):
            data = [((x - x0) / UPX, (y - ya) / (yb - ya)) for x, y in p]
            name = f"fig19{pn}_theta{th}_{kind}"
            with open(name + ".csv", "w") as f:
                w = csv.writer(f); w.writerow(["U_over_Uref", "r_over_d"]); [w.writerow([round(u, 4), round(r, 4)]) for u, r in data]
            out[name] = dict(points=len(p), r=[round(min(r for u, r in data), 3), round(max(r for u, r in data), 3)], Umax=round(max(u for u, r in data), 3))
            paths.append(p)
        pts = [((cx - x0) / UPX, (cy - ya) / (yb - ya)) for cx, cy in pk if x0 - 15 < cx < x1 + 40 - (j < 2) * 40 and ya - 5 < cy < yb + 5]
        with open(f"fig19{pn}_theta{th}_humphrey.csv", "w") as f:
            w = csv.writer(f); w.writerow(["U_over_Uref", "r_over_d"]); [w.writerow([round(u, 4), round(r, 4)]) for u, r in sorted(pts, key=lambda q: q[1])]
        out[f"fig19{pn}_theta{th}_humphrey"] = len(pts)
core = np.zeros_like(ink)
for a, b in PAN.values(): core[int(a) + 18:int(b) - 18, int(XL) + 18:int(XR) - 18] = True
cov, miss = coverage(curve, paths, core)
print(json.dumps(out, indent=0)); print(f"{len(pk)} points; curve ink explained {cov:.1%}")
im = Image.open(SRC).convert("RGB"); d = ImageDraw.Draw(im)
for k, p in enumerate(paths): d.line([tuple(q) for q in p], fill=((220, 0, 0) if k % 2 == 0 else (0, 110, 230)), width=2)
for cx, cy in pk: d.ellipse([cx - 6, cy - 6, cx + 6, cy + 6], outline=(0, 160, 160), width=3)
ys, xs = np.nonzero(miss); [d.point((x, y), fill=(255, 0, 255)) for x, y in zip(xs, ys)]
im.save("fig19_overlay.png")
