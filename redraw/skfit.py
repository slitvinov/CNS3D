import sys, os, subprocess, numpy as np
from PIL import Image
from scipy import ndimage as ndi, optimize
scan, dat, pdf = sys.argv[1:4]
env = dict(os.environ); POP = os.environ.get("PDFTOPPM", "pdftoppm")
A = 1 - np.asarray(Image.open(scan).convert("L"), np.float32) / 255
def z(m, s): b = ndi.gaussian_filter(m, s); return (b - b.mean()) / (b.std() + 1e-9)
def render(p):
    subprocess.run(["./sketchd", dat, pdf, f"{p[0]:.4f} {p[1]:.3f} {p[2]:.4f}"], env=env, capture_output=True)
    subprocess.run([POP, "-r", "300", "-gray", "-png", "-singlefile", pdf, pdf[:-4]], capture_output=True)
    r = 1 - np.asarray(Image.open(pdf[:-4] + ".png").convert("L"), np.float32) / 255
    o = np.zeros_like(A); h, w = min(r.shape[0], A.shape[0]), min(r.shape[1], A.shape[1]); o[:h, :w] = r[:h, :w]; return o
p = np.array([1.0, 0.0, 1.0])
for s, n in ((2.0, 60), (1.0, 60)):
    S = z(A, s); f = lambda q: -float((z(render(q), s) * S).mean())
    p = optimize.minimize(f, p, method="Nelder-Mead", options=dict(initial_simplex=[p, p + [0.15, 0, 0], p + [0, 3, 0], p + [0, 0, 0.2]], maxfev=n)).x
render(p); print(" ".join(f"{v:.4f}" for v in p))
lines = open(dat).read().split("\n"); nw, nh, nl, npo, nt = map(int, lines[0].split())
i = 1
for _ in range(nl): n = int(lines[i].split()[0]); i += 1 + n
for _ in range(npo): n = int(lines[i]); i += 1 + n
R = render(p); Ab = ndi.gaussian_filter(A, 1.0); Rb = ndi.gaussian_filter(R, 1.0); moved = 0
for _ in range(nt):
    x0, x1, yb, ch = map(float, lines[i].split()); y0, y1 = int(yb - ch * 1.3) - 8, int(yb + ch * 0.35) + 8; a0, a1 = int(x0) - 8, int(x1) + 8
    if y0 < 0 or a0 < 0 or y1 > A.shape[0] or a1 > A.shape[1]: i += 2; continue
    best = (-1e9, 0, 0); S_ = Ab[y0:y1, a0:a1]
    for dy in range(-6, 7):
        for dx in range(-6, 7):
            T_ = np.roll(np.roll(Rb, dy, 0), dx, 1)[y0:y1, a0:a1]; s_ = float((S_ * T_).sum())
            if s_ > best[0]: best = (s_, dx, dy)
    if best[1] or best[2]: moved += 1
    lines[i] = f"{x0 + best[1]:.1f} {x1 + best[1]:.1f} {yb + best[2]:.1f} {ch:.1f}"; i += 2
open(dat, "w").write("\n".join(lines)); render(p)
