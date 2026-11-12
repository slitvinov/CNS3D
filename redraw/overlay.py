import sys, json, numpy as np
from PIL import Image
from scipy import ndimage as ndi
scan, red, out = sys.argv[1:4]
A = np.asarray(Image.open(scan).convert("L")) < 128; B0 = np.asarray(Image.open(red).convert("L")) < 128
B = np.zeros_like(A); h, w = min(A.shape[0], B0.shape[0]), min(A.shape[1], B0.shape[1]); B[:h, :w] = B0[:h, :w]
da, db = ndi.distance_transform_edt(~A), ndi.distance_transform_edt(~B)
print(json.dumps(dict(scan_ink_within_2px=round(float(np.mean(db[A] <= 2)), 4), redraw_ink_within_2px=round(float(np.mean(da[B] <= 2)), 4))))
c = np.zeros(A.shape + (3,), np.uint8); c[A] = [255, 0, 255]; c[B] = [0, 255, 0]; c[A & B] = [255, 255, 255]
Image.fromarray(c).save(out + "_composite.png")
