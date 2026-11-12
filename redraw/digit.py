import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi

def load(path, thr=128):
    g = np.asarray(Image.open(path).convert("L"), np.float32)
    return g, g < thr

def runs(v, minlen):
    out, s = [], None
    for i, b in enumerate(list(v) + [False]):
        if b and s is None: s = i
        if not b and s is not None:
            if i - s >= minlen: out.append((s, i))
            s = None
    return out

def frame(ink, frac=0.5):
    H, W = ink.shape
    cols = [x for x in range(W) if any(e - s > frac * H for s, e in runs(ink[:, x], 1))]
    rows = [y for y in range(H) if any(e - s > frac * W for s, e in runs(ink[y, :], 1))]
    group = lambda v: [int(np.mean(g)) for g in np.split(np.array(v), np.where(np.diff(v) > 3)[0] + 1)] if v else []
    return group(cols), group(rows)

def ticks(ink, axis, at, lo, hi, depth=(3, 14), side=-1):
    d0, d1 = depth
    if axis == "x":
        band = ink[at + side * d1: at + side * d0, :] if side < 0 else ink[at + d0: at + d1, :]
        prof = band[:, lo:hi].mean(0)
    else:
        band = ink[:, at + d0: at + d1] if side > 0 else ink[:, at - d1: at - d0]
        prof = band[lo:hi, :].mean(1)
    on = prof > 0.6
    return [lo + (s + e - 1) / 2 for s, e in runs(on, 1) if e - s <= 8]

def calib(px, val):
    px, val = np.asarray(px, float), np.asarray(val, float)
    A = np.vstack([px, np.ones_like(px)]).T; (a, b), *_ = np.linalg.lstsq(A, val, rcond=None)
    res = (val - (a * px + b)) / a
    return a, b, float(np.sqrt(np.mean(res ** 2)))

def markers(gray, center, r, thr=0.75, mask=None):
    cx, cy = center; t = gray[cy - r: cy + r + 1, cx - r: cx + r + 1].astype(np.float32)
    t = (255 - t); t -= t.mean(); t /= np.linalg.norm(t) + 1e-9
    img = 255 - gray.astype(np.float32); k = t.shape[0]
    num = ndi.correlate(img, t, mode="constant")
    mean = ndi.uniform_filter(img, k, mode="constant"); sq = ndi.uniform_filter(img ** 2, k, mode="constant")
    den = np.sqrt(np.maximum(sq - mean ** 2, 1e-9) * k * k)
    cc = num / den
    if mask is not None: cc[~mask] = 0
    peaks = (cc == ndi.maximum_filter(cc, size=2 * r)) & (cc > thr)
    ys, xs = np.nonzero(peaks)
    return sorted(zip(xs.tolist(), ys.tolist(), cc[ys, xs].round(3).tolist()))

def refine_circle(ink, x, y, r):
    w = ink[y - r: y + r + 1, x - r: x + r + 1]; yy, xx = np.nonzero(w)
    return (x - r + xx.mean(), y - r + yy.mean()) if len(xx) else (x, y)

def trace(ink, way, step=3.0, maxgap=40, cone=0.6, pull=0.5):
    H, W = ink.shape; ys, xs = np.nonzero(ink); pts = np.stack([xs, ys], 1).astype(float)
    tree_idx = {}
    for i, (x, y) in enumerate(pts): tree_idx.setdefault((int(x) // 8, int(y) // 8), []).append(i)
    def near(p, rad):
        cx, cy = int(p[0]) // 8, int(p[1]) // 8; R = int(rad) // 8 + 1; idx = []
        for i in range(cx - R, cx + R + 1):
            for j in range(cy - R, cy + R + 1): idx += tree_idx.get((i, j), [])
        if not idx: return np.zeros((0, 2))
        q = pts[idx]; return q[np.hypot(*(q - p).T) <= rad]
    def centre(p, d):
        n = np.array([-d[1], d[0]]); q = near(p, 6)
        if not len(q): return p
        along = (q - p) @ d; acr = (q - p) @ n; sel = np.abs(along) <= 2
        return p + n * np.median(acr[sel]) if sel.any() else p
    path, found = [np.array(way[0], float)], [True]
    for target in way[1:]:
        target = np.array(target, float)
        d = target - path[-1]; d /= np.linalg.norm(d) + 1e-9
        while np.linalg.norm(target - path[-1]) > step * 1.5:
            p = path[-1]; tdir = target - p; tdir /= np.linalg.norm(tdir)
            want = (1 - pull) * d + pull * tdir; want /= np.linalg.norm(want)
            got = None
            for dist in np.arange(step, maxgap + step, step):
                q = near(p + want * dist, step)
                if len(q):
                    v = q - p; v /= np.linalg.norm(v, axis=1)[:, None] + 1e-9
                    ok = q[(v @ want) > np.cos(cone)]
                    if len(ok):
                        got = ok[np.argmin(np.hypot(*(ok - (p + want * dist)).T))]; break
            if got is None:
                nxt = p + want * step; f = False
            else:
                nxt = centre(got, want); f = dist <= step + 1e-6
                if dist > step:
                    for s in np.arange(step, dist, step): path.append(p + want * s); found.append(False)
            nd = nxt - path[-1]
            if np.linalg.norm(nd) > 1e-6: d = 0.6 * d + 0.4 * nd / np.linalg.norm(nd); d /= np.linalg.norm(d)
            path.append(nxt); found.append(f)
            if not (0 <= nxt[0] < W and 0 <= nxt[1] < H) or len(path) > 20000: break
    return np.array(path), np.array(found)

def ink_distance(ink, path):
    dt = ndi.distance_transform_edt(~ink)
    xs = np.clip(path[:, 0].round().astype(int), 0, ink.shape[1] - 1); ys = np.clip(path[:, 1].round().astype(int), 0, ink.shape[0] - 1)
    return dt[ys, xs]

def coverage(ink, paths, region, radius=2.5):
    m = np.zeros_like(ink)
    for p in paths:
        xs = np.clip(p[:, 0].round().astype(int), 0, ink.shape[1] - 1); ys = np.clip(p[:, 1].round().astype(int), 0, ink.shape[0] - 1)
        m[ys, xs] = True
    near = ndi.distance_transform_edt(~m) <= radius
    tgt = ink & region
    return float((tgt & near).sum() / max(tgt.sum(), 1)), tgt & ~near

def trace_fn(ink, guide, x_from, x_to, tol=6.0, blend=0.5, passes=2, smooth=25):
    gx, gy = np.array(guide, float).T
    for it in range(passes):
        G = lambda x: np.interp(x, gx, gy)
        xs, ys = [], []
        for x in range(int(x_from), int(x_to) + 1):
            col = ink[:, x]; rr = runs(col, 1)
            if not rr: continue
            cen = np.array([(s + e - 1) / 2 for s, e in rr]); wid = np.array([e - s for s, e in rr])
            pred = G(x)
            if len(xs) >= 4 and x - xs[-1] < 40:
                k = np.polyfit(xs[-6:], ys[-6:], 1); pred = blend * pred + (1 - blend) * np.polyval(k, x)
            j = np.argmin(np.abs(cen - pred))
            if abs(cen[j] - pred) <= tol + (wid[j] > 8) * 0 and wid[j] <= 14:
                xs.append(x); ys.append(cen[j])
        if passes > 1 and it < passes - 1 and len(xs) > smooth:
            ker = np.ones(smooth) / smooth; sy = np.convolve(ys, ker, "same")
            gx, gy = np.array(xs, float)[smooth // 2: -smooth // 2], sy[smooth // 2: -smooth // 2]
    return np.array(xs, float), np.array(ys, float)

def chain(mask, x0, x1, y0, y1, start, end=None, gap=60, cost=None, gapw=3.0):
    from scipy.spatial import cKDTree
    from skimage.morphology import skeletonize
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    sk = skeletonize(mask[int(y0):int(y1), int(x0):int(x1)]); ys, xs = np.nonzero(sk)
    P = np.stack([xs + int(x0), ys + int(y0)], 1).astype(float)
    if len(P) < 10: return np.zeros((0, 2))
    pr = np.array(list(cKDTree(P).query_pairs(gap)))
    if not len(pr): return np.zeros((0, 2))
    d = np.hypot(*(P[pr[:, 0]] - P[pr[:, 1]]).T); w = np.where(d <= 1.5, d, 1.5 + (d - 1.5) * gapw)
    if cost is not None:
        c = cost[P[:, 1].astype(int), P[:, 0].astype(int)]; w = w * (c[pr[:, 0]] + c[pr[:, 1]]) / 2
    G = coo_matrix((np.r_[w, w], (np.r_[pr[:, 0], pr[:, 1]], np.r_[pr[:, 1], pr[:, 0]])), shape=(len(P), len(P))).tocsr()
    i = int(np.argmin(np.hypot(*(P - start).T))); j = int(np.argmin(np.hypot(*(P - end).T)))
    dist, pred = dijkstra(G, indices=i, return_predecessors=True)
    if not np.isfinite(dist[j]): return np.zeros((0, 2))
    path = [j]
    while path[-1] != i: path.append(pred[path[-1]])
    return P[path[::-1]]
