import sys
from PIL import Image, ImageDraw
scan, red, comp, out = sys.argv[1:5]
ims = [Image.open(p).convert("RGB") for p in (scan, red, comp)]
h = max(i.height for i in ims); W = sum(i.width for i in ims) + 40
T = Image.new("RGB", (W, h + 50), "white"); x = 0; d = ImageDraw.Draw(T)
for im, lab in zip(ims, ("scan", "DISLIN redraw", "overlay: scan magenta, redraw green, both white")):
    T.paste(im, (x, 50)); d.text((x + 10, 15), lab, fill="black"); x += im.width + 20
T.thumbnail((3000, 3000)); T.save(out)
