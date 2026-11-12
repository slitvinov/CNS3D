import sys, re, json, subprocess, os, numpy as np, xml.etree.ElementTree as ET
from PIL import Image
from scipy import ndimage as ndi
from skimage.morphology import skeletonize
from skimage.measure import approximate_polygon, find_contours

TESS = os.environ.get("TESSERACT", "tesseract")
src, out = sys.argv[1], sys.argv[2]
fix = json.load(open(sys.argv[3])) if len(sys.argv) > 3 and os.path.exists(sys.argv[3]) else {}
boxes = fix if isinstance(fix, list) else None
if boxes is not None: fix = {}
g = np.asarray(Image.open(src).convert("L")); ink = g < 128; H, W = ink.shape

env = dict(os.environ, OMP_THREAD_LIMIT="1")
big = os.path.splitext(out)[0] + "_x3.png"
if boxes is not None: Image.new("L", (8, 8), 255).save(big); Image.open(src).convert("L").resize((3 * W, 3 * H), Image.LANCZOS).save(big)
def hocr(psm):
    h = subprocess.run([TESS, big, "-", "--psm", str(psm), "hocr"], capture_output=True, text=True, env=env).stdout
    return ET.fromstring(h)
def bb(e): return [int(v) for v in re.search(r"bbox (\d+) (\d+) (\d+) (\d+)", e.get("title")).groups()]
def ok(t): return re.search(r"[A-Za-z0-9]", t)
words = []
def overlaps(b): return any(b[0] < w[2] and b[2] > w[0] and b[1] < w[3] and b[3] > w[1] for w in words)
for ln in hocr(3).iter():
    if ln.get("class") not in ("ocr_line", "ocr_caption", "ocr_textfloat", "ocr_header"): continue
    ws = []
    for e in ln.iter():
        if e.get("class") != "ocrx_word": continue
        t = "".join("".join(e.itertext()).split()); c = int(re.search(r"x_wconf (\d+)", e.get("title")).group(1))
        if t and ok(t) and c >= fix.get("conf", 40): ws.append(bb(e) + [t])
    grp = []
    for w in ws:
        if grp and w[0] - grp[-1][-1][2] < 1.6 * (w[3] - w[1]): grp[-1].append(w)
        else: grp.append([w])
    for g_ in grp:
        x0 = min(w[0] for w in g_) // 3; y0 = min(w[1] for w in g_) // 3; x1 = (max(w[2] for w in g_) + 2) // 3; y1 = (max(w[3] for w in g_) + 2) // 3
        if 7 <= y1 - y0 <= 70 and x1 - x0 < 0.6 * W: words.append([x0, y0, x1, y1, " ".join(w[4] for w in g_)])
for e in hocr(11).iter():
    if e.get("class") != "ocrx_word": continue
    t = "".join("".join(e.itertext()).split()); c = int(re.search(r"x_wconf (\d+)", e.get("title")).group(1))
    x0, y0, x1, y1 = bb(e); b = [x0 // 3, y0 // 3, (x1 + 2) // 3, (y1 + 2) // 3]
    if t and ok(t) and c >= 75 and 7 <= b[3] - b[1] <= 70 and not overlaps(b): words.append(b + [t])
for b in fix.get("drop", []):
    words = [w for w in words if not (w[0] < b[2] and w[2] > b[0] and w[1] < b[3] and w[3] > b[1])]
for w in words:
    if w[4] in fix.get("text", {}): w[4] = fix["text"][w[4]]
words += fix.get("add", [])
if boxes is not None: words = [list(b) for b in boxes]
words = [w if len(w) == 7 else list(w[:5]) + [w[3], w[3] - w[1]] for w in words]
tmask = np.zeros_like(ink)
for x0, y0, x1, y1, t, *_ in words: tmask[max(y0 - 2, 0):y1 + 3, max(x0 - 2, 0):x1 + 3] = True
art = ink & ~tmask
dt = ndi.distance_transform_edt(art)
sk0 = skeletonize(art); hw = float(np.median(dt[sk0])) if sk0.any() else 1.5
r = max(2, int(round(hw + 1.6)))
solid = ndi.binary_opening(art, structure=np.ones((3, 3)), iterations=r)
solid = ndi.binary_dilation(solid, iterations=1) & art
lab, n = ndi.label(solid); sz = ndi.sum(solid, lab, range(1, n + 1))
polys = []
for k in range(1, n + 1):
    if sz[k - 1] < 25: solid[lab == k] = False; continue
    m = np.pad(lab == k, 1)
    for c in find_contours(m.astype(float), 0.5):
        c = approximate_polygon(c, 0.7)
        if len(c) >= 3: polys.append([(float(x - 1), float(y - 1)) for y, x in c])
thin = art & ~solid
sk = skeletonize(thin)
nb = ndi.convolve(sk.astype(int), np.ones((3, 3), int), mode="constant") - 1; nb[~sk] = 0
node = sk & (nb != 2)
seen = np.zeros_like(sk); lines = []
off = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
def walk(y, x, py, px):
    p = [(y, x)]
    while True:
        seen[y, x] = True
        nx = [(y + a, x + b) for a, b in off if 0 <= y + a < H and 0 <= x + b < W and sk[y + a, x + b] and (y + a, x + b) != (py, px) and not (seen[y + a, x + b] and not node[y + a, x + b])]
        nx = [q for q in nx if q not in p[-3:]]
        if not nx: return p
        q = min(nx, key=lambda q: abs(q[0] - y) + abs(q[1] - x))
        p.append(q); py, px = y, x; y, x = q
        if node[y, x]: return p
for y, x in zip(*np.nonzero(node)):
    for a, b in off:
        yy, xx = y + a, x + b
        if 0 <= yy < H and 0 <= xx < W and sk[yy, xx] and not seen[yy, xx]:
            seen[y, x] = True; lines.append([(y, x)] + walk(yy, xx, y, x))
for y, x in zip(*np.nonzero(sk & ~seen)):
    if not seen[y, x]: lines.append(walk(y, x, -9, -9))
ws = []; L2 = []
for p in lines:
    if len(p) < 3: continue
    a = np.array(p, float); w = 2 * float(np.median(dt[a[:, 0].astype(int), a[:, 1].astype(int)]))
    a = approximate_polygon(a, 0.8); L2.append([(float(x), float(y)) for y, x in a]); ws.append(max(w, 1.0))
with open(out, "w") as f:
    f.write(f"{W} {H} {len(L2)} {len(polys)} {len(words)}\n")
    for p, w in zip(L2, ws): f.write(f"{len(p)} {w:.2f}\n" + "".join(f"{x:.1f} {y:.1f}\n" for x, y in p))
    for p in polys: f.write(f"{len(p)}\n" + "".join(f"{x:.1f} {y:.1f}\n" for x, y in p))
    for x0, y0, x1, y1, t, yb, ch in words: f.write(f"{x0} {x1} {yb} {ch}\n{t}\n")
print(json.dumps(dict(size=[W, H], lines=len(L2), filled=len(polys), labels=len(words), stroke=round(2 * hw, 2),
                      text=[w[4] for w in words])))
