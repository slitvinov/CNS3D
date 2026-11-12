import sys, json, csv; sys.path.insert(0, ".")
from digit import *
from scipy.signal import find_peaks
def lines(ink, axis, n=31):
    s = ndi.uniform_filter1d(ink.sum(axis).astype(float), 3)
    p, _ = find_peaks(s, distance=12, prominence=0.1 * s.max())
    return np.array(sorted(p, key=lambda i: -s[i])[:n] if len(p) > n else p, float)
def ends(ink, along, at, width=24, n=31):
    a = int(at); strip = ink[:, max(a - width // 2, 0):a + width // 2] if along == 0 else ink[max(a - width // 2, 0):a + width // 2, :]
    s = ndi.uniform_filter1d(strip.sum(1 - along).astype(float), 3)
    p, _ = find_peaks(s, distance=10, prominence=0.15 * s.max())
    p = np.sort(sorted(p, key=lambda i: -s[i])[:n]).astype(float)
    return p
PAD = 12
def pair(A, B, n=31, tol=8):
    swap = len(A) == n and len(B) != n
    R, O = (A, B) if swap else (B, A)
    assert len(R) == n, (len(A), len(B))
    d = np.median([O[np.argmin(np.abs(O - r))] - r for r in R])
    Q = np.array([O[np.argmin(np.abs(O - r - d))] if np.min(np.abs(O - r - d)) <= tol else r + d for r in R])
    return (R, Q) if swap else (Q, R)
def mesh(src):
    g, ink = load(src); ink = np.pad(ink, PAD); H, W = ink.shape
    vx = np.sort(lines(ink, 0)); assert len(vx) == 31, len(vx)
    VL = []
    for i in range(31):
        ys, xs = trace_fn(ink.T, [(0, vx[i]), (H - 1, vx[i])], PAD + 2, H - PAD - 2, tol=4, blend=0.8, passes=2, smooth=9)
        VL.append((np.polyfit(ys, xs, 2), ys, xs))
    m = ink.copy()
    for k, _, _ in VL:
        Y = np.arange(H); x = np.polyval(k, Y)
        for d in range(-3, 4): m[Y, np.clip((x + d).round().astype(int), 0, W - 1)] = False
    C = np.full((31, 31), np.nan); Y = np.arange(H)
    for i, (k, yv, _) in enumerate(VL):
        x = np.polyval(k, Y); prof = np.zeros(H)
        for d in list(range(-12, -4)) + list(range(5, 13)):
            prof += m[Y, np.clip((x + d).round().astype(int), 0, W - 1)]
        prof = ndi.uniform_filter1d(prof, 3); p, _ = find_peaks(prof, distance=6, prominence=0.2 * prof.max())
        p = p[(p > max(PAD, yv.min() - 4)) & (p < min(H - PAD, yv.max() + 4))]
        if len(p) > 31: p = np.sort(p[np.argsort(-prof[p])[:31]])
        if len(p) == 31 and abs(p[-1] - p[0]) > 0: C[i] = p
    ok = np.isfinite(C[:, 0]); ii = np.arange(31)
    def robust(v):
        u = ok.copy()
        for _ in range(3):
            k = np.polyfit(ii[u], v[u], 1); u = ok & (np.abs(v - np.polyval(k, ii)) < 6)
        return np.polyval(k, ii)
    top = robust(C[:, 0]); bot = robust(C[:, -1]); C[:] = np.nan
    for i, (k, yv, _) in enumerate(VL):
        x = np.polyval(k, Y); prof = np.zeros(H)
        for d in list(range(-12, -4)) + list(range(5, 13)):
            prof += m[Y, np.clip((x + d).round().astype(int), 0, W - 1)]
        prof = ndi.uniform_filter1d(prof, 3); p, _ = find_peaks(prof, distance=3, prominence=0.12 * prof.max())
        p = p[(p > top[i] - 5) & (p < bot[i] + 5)]
        if len(p) > 31: p = np.sort(p[np.argsort(-prof[p])[:31]])
        if len(p) == 31: C[i] = p
    good = np.isfinite(C[:, 0])
    gi = np.nonzero(good)[0]; bi = np.nonzero(~good)[0]
    for j in range(31):
        C[bi, j] = np.interp(bi, gi, C[gi, j])
        for i in bi[bi < gi[0]]: C[i, j] = C[gi[0], j] + (C[gi[0], j] - C[gi[1], j]) * (gi[0] - i) / (gi[1] - gi[0])
        for i in bi[bi > gi[-1]]: C[i, j] = C[gi[-1], j] + (C[gi[-1], j] - C[gi[-2], j]) * (i - gi[-1]) / (gi[-1] - gi[-2])
    for i in bi:
        k = VL[i][0]; x = np.polyval(k, Y); prof = np.zeros(H)
        for d in list(range(-12, -4)) + list(range(5, 13)):
            prof += m[Y, np.clip((x + d).round().astype(int), 0, W - 1)]
        prof = ndi.uniform_filter1d(prof, 3); e = C[i].copy()
        for j in range(31):
            h = 0.45 * min(e[j] - e[j - 1] if j else 99, e[j + 1] - e[j] if j < 30 else 99, 12)
            a_, b_ = int(round(e[j] - h)), int(round(e[j] + h)) + 1
            if b_ > a_ and prof[a_:b_].max() > 0.3 * prof.max(): C[i, j] = a_ + np.argmax(prof[a_:b_])
    st0 = dict(estimated_columns=bi.tolist())
    HL = []
    for j in range(31):
        gx = [np.polyval(VL[i][0], C[i, j]) for i in range(31)]; gp = list(zip(gx, C[:, j]))
        xs, ys = trace_fn(m, gp, gx[0] + 1, gx[-1] - 1, tol=3, blend=0.9, passes=1)
        k = np.polyfit(xs, ys, 3) if len(xs) > 40 else np.polyfit(gx, C[:, j], 3)
        HL.append((k, xs, ys))
    N = np.zeros((31, 31, 2))
    for i, (kv, _, _) in enumerate(VL):
        N[i, :, 1] = C[i] - PAD; N[i, :, 0] = np.polyval(kv, C[i]) - PAD
    for j, (k, xs, ys) in enumerate(HL): HL[j] = (None, xs, ys - np.interp(xs, N[:, j, 0] + PAD, N[:, j, 1] + PAD))
    res_h = np.concatenate([ys for k, xs, ys in HL]); res_v = np.concatenate([xs - np.polyval(k, ys) for k, ys, xs in VL])
    return ink, N, HL, VL, dict(rms_across=round(float(np.sqrt(np.mean(res_h ** 2))), 2), rms_up=round(float(np.sqrt(np.mean(res_v ** 2))), 2),
                               samples=int(len(res_h) + len(res_v)), **st0)
out = {}
for tag in ("a", "b"):
    f = f"fig5{tag}_page.png"
    ink, N, HL, VL, st = mesh(f)
    ang = float(np.mean([np.polyval(np.polyder(k), np.mean(ys)) for k, ys, xs in VL])); c0 = N.reshape(-1, 2).mean(0)
    R = np.array([[1, -ang], [ang, 1]]) / np.sqrt(1 + ang ** 2); N = (N - c0) @ R.T + c0; st["skew_deg"] = round(float(np.degrees(np.arctan(ang))), 3)
    x0, x1 = N[0, :, 0].mean(), N[-1, :, 0].mean(); y0, y1 = N[:, -1, 1].mean(), N[:, 0, 1].mean()
    with open(f"fig5{tag}_nodes.csv", "w") as fh:
        w = csv.writer(fh); w.writerow(["i", "j", "x_px", "y_px", "x", "y"])
        for i in range(31):
            for j in range(31):
                w.writerow([i + 1, j + 1, round(N[i, j, 0], 2), round(N[i, j, 1], 2), round((N[i, j, 0] - x0) / (x1 - x0), 5), round((y0 - N[i, j, 1]) / (y0 - y1), 5)])
    dx = np.diff(N[:, :, 0], axis=0); dy = -np.diff(N[:, :, 1], axis=1)
    st.update(box=[round(v, 1) for v in (x0, y1, x1, y0)], dx_range=[round(float(dx.min()), 1), round(float(dx.max()), 1)],
              dy_range=[round(float(dy.min()), 1), round(float(dy.max()), 1)])
    out[tag] = st
    N = (N - c0) @ R + c0
    im = Image.open(f).convert("RGB"); d = ImageDraw.Draw(im)
    for i in range(31): d.line([tuple(p) for p in N[i]], fill=(0, 110, 230), width=1)
    for j in range(31): d.line([tuple(p) for p in N[:, j]], fill=(220, 0, 0), width=1)
    im.save(f"fig5{tag}_overlay.png")
print(json.dumps(out, indent=0))
