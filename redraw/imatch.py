import sys, json, numpy as np
from PIL import Image
from scipy import ndimage as ndi, optimize
from skimage.registration import phase_cross_correlation
from skimage.metrics import structural_similarity as ssim
scan, red, out = sys.argv[1:4]
A = 1 - np.asarray(Image.open(scan).convert("L"), np.float32) / 255
B0 = 1 - np.asarray(Image.open(red).convert("L"), np.float32) / 255
H, W = max(A.shape[0], B0.shape[0]), max(A.shape[1], B0.shape[1])
pad = lambda m: np.pad(m, ((0, H - m.shape[0]), (0, W - m.shape[1])))
A, B0 = pad(A), pad(B0)
blur = lambda m: ndi.gaussian_filter(m, 2.0)
def ncc(a, b):
    a = a - a.mean(); b = b - b.mean(); return float((a * b).sum() / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))
def warp(m, p):
    sx, sy, sh, th, tx, ty = p
    c, s = np.cos(th), np.sin(th); M = np.array([[sy * c, -s + sh], [s, sx * c]])
    ctr = np.array([H / 2, W / 2]); off = ctr - M @ ctr + np.array([ty, tx])
    return ndi.affine_transform(m, M, offset=off, order=1)
def measures(a, b):
    ia, ib = a > 0.5, b > 0.5
    da, db = ndi.distance_transform_edt(~ia), ndi.distance_transform_edt(~ib)
    return dict(ncc=round(ncc(a, b), 4), ncc_blur2=round(ncc(blur(a), blur(b)), 4), ssim=round(float(ssim(1 - a, 1 - b, data_range=1.0)), 4),
                scan_ink_within_2px=round(float(np.mean(db[ia] <= 2)), 4), redraw_ink_within_2px=round(float(np.mean(da[ib] <= 2)), 4))
before = measures(A, B0)
shift, _, _ = phase_cross_correlation(blur(A), blur(B0), upsample_factor=10)
p0 = np.array([1, 1, 0, 0, shift[1], shift[0]], float)
f = lambda p: -ncc(blur(A), blur(warp(B0, p)))
res = optimize.minimize(lambda p: f(np.r_[p[:4], p[4:]]), np.array([1, 1, 0, 0, -shift[1], -shift[0]]), method="Powell",
                        options=dict(xtol=1e-4, ftol=1e-6, maxfev=4000))
p = res.x; B = warp(B0, p)
after = measures(A, B)
print(json.dumps(dict(before=before, after=after, registration=dict(scale_x=round(p[0], 4), scale_y=round(p[1], 4), shear=round(p[2], 5),
      rotation_deg=round(float(np.degrees(p[3])), 4), shift_x_px=round(p[4], 2), shift_y_px=round(p[5], 2))), indent=1))
ia, ib = A > 0.5, B > 0.5
comp = np.zeros((H, W, 3), np.uint8); comp[ia] = [255, 0, 255]; comp[ib] = [0, 255, 0]; comp[ia & ib] = [255, 255, 255]
Image.fromarray(comp).save(out + "_composite.png")
Image.fromarray((255 - np.abs(A - B) * 255).clip(0, 255).astype(np.uint8)).save(out + "_difference.png")
Image.fromarray(((1 - B) * 255).clip(0, 255).astype(np.uint8)).save(out + "_registered.png")
