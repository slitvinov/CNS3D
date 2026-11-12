import sys, json, csv; sys.path.insert(0, ".")
from digit import *
from scipy.signal import find_peaks
SRC = "fig16_deskew.png"
g, ink = load(SRC); H, W = ink.shape
cx, cy, ri, ro, q = np.load("fig16_ellipse.npy"); d = ro - ri
def peaks(s, dist, pr): s = ndi.uniform_filter1d(s.astype(float), 3); p, _ = find_peaks(s, distance=dist, prominence=pr * s.max()); return p
E = []
for x0, x1 in ((380, 620), (650, 860), (890, 1050)):
    r = int(cy + ri) - 15 + peaks(ink[int(cy + ri) - 15:int(cy + ro) + 15, x0:x1].sum(1), 4, 0.15); E.append((r - r[0]) / (r[-1] - r[0]))
for deg in (25, 45, 65, 85):
    rr = np.arange(ri - 15, ro + 15, 0.5); v = 0
    for off in (-6, 6):
        t = np.radians(deg) + off / ((ri + ro) / 2); v = v + ndi.map_coordinates(ink.astype(float), [cy + rr * np.cos(t), cx + q * rr * np.sin(t)], order=1)
    r = rr[peaks(v, 6, 0.15)]; E.append((r - r[0]) / (r[-1] - r[0]))
E = np.array([e for e in E if len(e) == 17]); eta = E.mean(0); eta_sd = E.std(0)
ins = 330 + peaks(ink[int(cy + ri) + 15:int(cy + ro) - 15, 330:1240].sum(0), 20, 0.4)
outs = 40 + peaks(ink[40:int(cy) + 60, int(cx + q * ri) + 15:int(cx + q * ro) - 15].sum(1), 20, 0.4)
s_in = np.sort((cx - ins) / (q * d))[::-1]; s_in = s_in[s_in > 0.05]
s_out = np.sort((cy - outs) / d); s_out = s_out[s_out > 0.05]
th = np.arange(0, 91, 10.0)
xs_ = json.load(open("fig16_xsec_px.json")); zc = np.array(xs_["z"], float); rc = np.array(xs_["r"], float)
zeta = (zc - zc[0]) / (zc[-1] - zc[0]); eta_x = (rc - rc[0]) / (rc[-1] - rc[0])
geo = dict(centre_px=[round(cx, 2), round(cy, 2)], Ri_px=round(ri, 2), d_px=round(d, 2), x_scale=round(q, 4), Ri_over_d=round(ri / d, 4),
           inlet_stations_d=np.round(s_in, 4).tolist(), bend_stations_deg=th.tolist(), outlet_stations_d=np.round(s_out, 4).tolist(),
           eta=np.round(eta, 4).tolist(), eta_spread=round(float(eta_sd.max()), 4), xsec_zeta=np.round(zeta, 4).tolist(), xsec_eta=np.round(eta_x, 4).tolist(),
           xsec_box_px=[float(zc[0]), float(rc[0]), float(zc[-1]), float(rc[-1])], xsec_width_over_d=round(float((zc[-1] - zc[0]) / (q * d)), 4))
json.dump(geo, open("fig16_grid.json", "w"), indent=0)
with open("fig16_nodes.csv", "w") as f:
    w = csv.writer(f); w.writerow(["i", "k", "x_over_d", "y_over_d"]); i = 0
    st = [("in", s) for s in s_in] + [("bend", t) for t in th] + [("out", s) for s in s_out]
    for kind, v in st:
        i += 1
        for k, e in enumerate(eta):
            r = ri / d + e
            if kind == "in": x, y = -v, -r
            elif kind == "bend": x, y = r * np.sin(np.radians(v)), -r * np.cos(np.radians(v))
            else: x, y = r, v
            w.writerow([i, k + 1, round(x, 5), round(y, 5)])
print(json.dumps(geo)); print("stations", len(s_in) + len(th) + len(s_out), "x", len(eta))
