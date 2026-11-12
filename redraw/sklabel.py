import sys, json, numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from skimage.feature import match_template

F = {"": "/usr/share/fonts/urw-base35/NimbusSans-Bold.otf", "@": "/usr/share/fonts/urw-base35/StandardSymbolsPS.t1"}
GR = dict(zip("zxhepqrDnmsa", "ζξηεπθρΔνμσα"))
src, spec, out = sys.argv[1:4]
g = 1 - np.asarray(Image.open(src).convert("L"), np.float32) / 255; H, W = g.shape
gb = ndi.gaussian_filter(g, 0.8)
def render(t, px, wf):
    font = ImageFont.truetype(F["@"] if t[0] == "@" else F[""], px); s = "".join(GR.get(c, c) for c in t[1:]) if t[0] == "@" else t
    l, tp, r, b = font.getbbox(s, anchor="ls"); im = Image.new("L", (r - l + 8, b - tp + 8), 0)
    ImageDraw.Draw(im).text((4 - l, 4 - tp), s, fill=255, font=font, anchor="ls")
    im = im.resize((max(3, int(round(im.width * wf))), im.height), Image.LANCZOS)
    a = np.asarray(im, np.float32) / 255; ys, xs = np.nonzero(a > 0.5)
    cap = -ImageFont.truetype(F[""], px).getbbox("H", anchor="ls")[1]
    render.base = 4 - tp - ys.min(); render.cap = cap
    return a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
res = []
for t, cx, cy, h in json.load(open(spec)):
    best = (-9, None)
    for hh in np.arange(0.8, 1.26, 0.05) * h:
        for wf in (0.6, 0.68, 0.76, 0.84, 0.92, 1.0):
            T = render(t, max(6, int(round(hh * 1.38))), wf)
            th, tw = T.shape; R = 22
            x0, y0 = int(cx - tw / 2 - R), int(cy - th / 2 - R)
            x0c, y0c = max(x0, 0), max(y0, 0); x1c, y1c = min(x0 + tw + 2 * R, W), min(y0 + th + 2 * R, H)
            win = gb[y0c:y1c, x0c:x1c]
            if win.shape[0] <= th or win.shape[1] <= tw: continue
            m = match_template(win, ndi.gaussian_filter(T, 0.8)); m = np.where(np.isfinite(m) & (m <= 1.0), m, -1)
            from scipy.signal import fftconvolve
            mass = fftconvolve(win, np.ones_like(T), mode="valid"); m = np.where((mass > 0.6 * T.sum()) & (mass < 1.8 * T.sum()), m, -1)
            k = np.unravel_index(m.argmax(), m.shape); s = m[k]
            if s > best[0]: best = (s, [x0c + k[1], y0c + k[0], x0c + k[1] + tw, y0c + k[0] + th, y0c + k[0] + render.base, render.cap])
    b = best[1]; res.append([int(v) for v in b[:4]] + [t, round(float(b[4]), 1), round(float(b[5]), 1)]); print(f"{t:28s} {best[0]:.2f} {[round(float(v), 1) for v in best[1]]}")
json.dump(res, open(out, "w"))
