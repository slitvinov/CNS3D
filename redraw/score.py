import sys, os, json, subprocess, re
import numpy as np
from PIL import Image
from scipy import ndimage as ndi, optimize
POP = os.environ.get("PDFTOPPM", "pdftoppm")
FIGS = [("1", "../fig/p16-1.png", "fig1_dislin.pdf"), ("2", "../fig/p18-1.png", "fig2_dislin.pdf"), ("3", "../fig/p21-1.png", "fig3_dislin.pdf"),
        ("4", "../fig/p22-1.png", "fig4_dislin.pdf"), ("5a", "fig5a_deskew.png", "fig5a_dislin.pdf"), ("5b", "fig5b_deskew.png", "fig5b_dislin.pdf"),
        ("6", "fig6_page.png", "fig6_dislin.pdf"), ("7", "fig7_deskew.png", "fig7_dislin.pdf"), ("8", "fig8_deskew.png", "fig8_dislin.pdf"),
        ("9", "../fig/p26-1.png", "fig9_dislin.pdf"), ("10", "../fig/p27-1.png", "fig10_dislin.pdf"), ("11", "fig11_deskew.png", "fig11_dislin.pdf"),
        ("12", "../fig/p28-1.png", "fig12_dislin.pdf"), ("13", "fig13_deskew.png", "fig13_dislin.pdf"), ("14", "fig14_deskew.png", "fig14_dislin.pdf"),
        ("15", "fig15_deskew.png", "fig15_dislin.pdf"), ("16a", "fig16a_scan.png", "fig16a_dislin.pdf"), ("16b", "fig16b_deskew.png", "fig16_dislin.pdf"),
        ("17", "fig17_page.png", "fig17_dislin.pdf"), ("18", "fig18_page.png", "fig18_dislin.pdf"), ("19", "fig19_deskew.png", "fig19_dislin.pdf"),
        ("A-1", "../fig/p38-1.png", "figA1_dislin.pdf"), ("A-2", "../fig/p39-1.png", "figA2_dislin.pdf")]
CROP = {"16b": (1000, None)}
TILT = {"14": ("../fig/p29-2.png", 0.311), "15": ("../fig/p29-3.png", 0.287)}
def ink(im): return np.asarray(im.convert("L")) < 128
def prep(m, f): return ndi.gaussian_filter(m[::f, ::f].astype(np.float32), 1.5)
def warp(b, t, dx, dy):
    c = np.array(b.shape) / 2.0; a = np.radians(t); R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    off = c - R @ c - np.array([dy, dx]); return ndi.affine_transform(b, R, offset=off, order=1)
def ncc(a, b): a = a - a.mean(); b = b - b.mean(); return float((a * b).sum() / (np.sqrt((a * a).sum() * (b * b).sum()) + 1e-12))
def rigid(A, B, f=2):
    a, b = prep(A, f), prep(B, f)
    best = None
    for t0 in (-0.6, -0.3, 0.0, 0.3, 0.6):
        r = optimize.minimize(lambda p: -ncc(a, warp(b, p[0], p[1], p[2])), [t0, 0, 0], method="Powell", options=dict(xtol=1e-3, ftol=1e-7))
        if best is None or r.fun < best.fun: best = r
    t, dx, dy = best.x; return float(t), float(dx * f), float(dy * f), -float(best.fun)
def pair(scan, pdf, tag):
    out = f"/tmp/score_{tag}"
    subprocess.run([POP, "-r", "300", "-gray", "-png", "-singlefile", pdf, out], check=True)
    A = ink(Image.open(scan)); R = ink(Image.open(out + ".png"))
    return A, R
def measures(A, B):
    dA, dB = ndi.distance_transform_edt(~A), ndi.distance_transform_edt(~B)
    sa, rb = dB[A], dA[B]
    return dict(scan2=float(np.mean(sa <= 2)), redraw2=float(np.mean(rb <= 2)), scan1=float(np.mean(sa <= 1)), redraw1=float(np.mean(rb <= 1)),
                scan_mean_px=float(sa.mean()), redraw_mean_px=float(rb.mean()), scan_p95_px=float(np.percentile(sa, 95)), redraw_p95_px=float(np.percentile(rb, 95)))
def score(n, scan, pdf):
    A, R = pair(scan, pdf, n)
    flags = []
    if abs(A.shape[0] - R.shape[0]) > 2 or abs(A.shape[1] - R.shape[1]) > 2: flags.append(f"size scan {A.shape[1]}x{A.shape[0]} redraw {R.shape[1]}x{R.shape[0]}")
    B = np.zeros_like(A); h, w = min(A.shape[0], R.shape[0]), min(A.shape[1], R.shape[1]); B[:h, :w] = R[:h, :w]
    if n in CROP: A, B = A[CROP[n][0]:CROP[n][1]], B[CROP[n][0]:CROP[n][1]]
    m = measures(A, B)
    m["floor_scan2"] = measures(A, np.roll(B, (10, 10), (0, 1)))["scan2"]
    t, dx, dy, q = rigid(A, B); m["rot_rel_deg"] = t; m["shift_px"] = [dx, dy]
    if abs(t) > 0.1: flags.append(f"redraw must turn {t:+.2f} deg to fit the scan")
    if np.hypot(dx, dy) > 1.5: flags.append(f"redraw must shift {dx:+.1f},{dy:+.1f} px to fit the scan")
    if m["scan2"] - m["floor_scan2"] < 0.3: flags.append("score close to the chance floor")
    m["flags"] = flags
    overlay(A, B, f"compare/fig{n}.png")
    if n in TILT:
        sc, ang = TILT[n]; A2 = ink(Image.open(sc)); R2 = Image.fromarray(np.where(B, 0, 255).astype(np.uint8)).rotate(ang, resample=Image.BICUBIC, fillcolor=255)
        B2 = np.zeros_like(A2); r2 = ink(R2); h, w = min(A2.shape[0], r2.shape[0]), min(A2.shape[1], r2.shape[1]); B2[:h, :w] = r2[:h, :w]
        t = measures(A2, B2); m["tilted"] = dict(scan2=t["scan2"], redraw2=t["redraw2"]); m["level_on_printed"] = {k: measures(A2, np.pad(B, ((0, max(0, A2.shape[0] - B.shape[0])), (0, max(0, A2.shape[1] - B.shape[1]))))[:A2.shape[0], :A2.shape[1]])[k] for k in ("scan2", "redraw2")}
        overlay(A2, B2, f"compare/fig{n}_tilt.png")
    return m
def overlay(A, B, out):
    o = np.full(A.shape + (3,), 255, np.uint8); o[A] = (160, 160, 160); o[B] = (215, 25, 28); o[A & B] = (0, 0, 0)
    S = Image.fromarray(o); S.thumbnail((1400, 1400), Image.LANCZOS); S.quantize(colors=64, dither=Image.Dither.NONE).save(out, optimize=True)
def pct(v): return f"{100 * v:.0f}"
def check(res):
    bad = []
    g = open("GALLERY.md").read()
    for n, m in res.items():
        sec = re.search(r"### Figure " + re.escape(n) + r"\..*?\nMatch: ([0-9.]+)(?: / ([0-9.]+))?", g, re.S)
        if not sec: bad.append(f"{n}: no Match line in GALLERY.md"); continue
        a, b = float(sec.group(1)), sec.group(2)
        if abs(a - 100 * m["scan2"]) > 0.6 or (b and abs(float(b) - 100 * m["redraw2"]) > 0.6): bad.append(f"{n}: GALLERY says {sec.group(1)} / {b}, measured {pct(m['scan2'])} / {pct(m['redraw2'])}")
    return bad
if __name__ == "__main__":
    only = sys.argv[1:]
    res = json.load(open("scores.json")) if os.path.exists("scores.json") else {}
    if only == ["check"]:
        bad = check(res); print("\n".join(bad) if bad else "GALLERY.md agrees with scores.json"); sys.exit(1 if bad else 0)
    for n, scan, pdf in FIGS:
        if only and n not in only: continue
        m = score(n, scan, pdf); res[n] = m
        print(f"{n:4s} {pct(m['scan2']):>3s} / {pct(m['redraw2']):>3s}  1px {pct(m['scan1']):>3s} / {pct(m['redraw1']):>3s}  mean {m['scan_mean_px']:.2f} / {m['redraw_mean_px']:.2f} px  p95 {m['scan_p95_px']:.1f} / {m['redraw_p95_px']:.1f}  floor {pct(m['floor_scan2'])}  rot {m['rot_rel_deg']:+.2f}  {'; '.join(m['flags'])}", flush=True)
    json.dump(res, open("scores.json", "w"), indent=1, sort_keys=True)
