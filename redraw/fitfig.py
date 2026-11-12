import os, sys, json, subprocess, numpy as np
from PIL import Image
from scipy import ndimage as ndi, optimize
F = sys.argv[1]; cfg = json.load(open(f"fig{F}_fit.json"))
env = dict(os.environ)
POP = os.environ.get("PDFTOPPM", "pdftoppm")
scan = 1 - np.asarray(Image.open(cfg["scan"]).convert("L"), np.float32) / 255
names = cfg["names"]; p = dict(zip(names, [float(v) for v in open(f"fig{F}_style.dat").read().split()]))
fit = cfg.get("fit", names); step = cfg["step"]
def render(q):
    pp = dict(p); pp.update(zip(fit, q))
    open(f"fig{F}_style.dat", "w").write(" ".join(f"{pp[n]:.4f}" for n in names) + "\n")
    subprocess.run([f"./fig{F}d"], env=env, capture_output=True)
    subprocess.run([POP, "-r", "300", "-gray", "-png", "-singlefile", f"fig{F}_dislin.pdf", f"fig{F}_dislin"], capture_output=True)
    r = 1 - np.asarray(Image.open(f"fig{F}_dislin.png").convert("L"), np.float32) / 255
    o = np.zeros_like(scan); h, w = min(r.shape[0], scan.shape[0]), min(r.shape[1], scan.shape[1]); o[:h, :w] = r[:h, :w]; return o
def z(m, s): b = ndi.gaussian_filter(m, s); return (b - b.mean()) / (b.std() + 1e-9)
q = np.array([p[n] for n in fit]); st = np.array([step[n] for n in fit]); n = 0
for sig, nfev in cfg.get("stages", [[6, 400], [3, 400], [1.5, 400]]):
    S = z(scan, sig); best = [-9, None]
    def f(x):
        global n; n += 1; s = float((z(render(x), sig) * S).mean())
        if s > best[0]: best[0], best[1] = s, x.copy()
        return -s
    optimize.minimize(f, q, method="Nelder-Mead", options=dict(initial_simplex=np.vstack([q] + [q + np.eye(len(q))[i] * st[i] for i in range(len(q))]), maxfev=nfev, xatol=0.02, fatol=1e-6))
    q = best[1]; st = st / 2; print(f"blur {sig}: NCC {best[0]:.4f}", flush=True)
p.update(zip(fit, q)); render(q)
print(f"{n} renders"); print(json.dumps({k: round(v, 3) for k, v in p.items()}))
