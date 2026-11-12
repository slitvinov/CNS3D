import csv, numpy as np, sys
n = 31; NZ = 4; dz = 1.0 / (n - 1); f = sys.argv[1]
X = np.zeros((n, n)); Y = np.zeros((n, n))
if f.startswith("warp"):
    a_ = float(f[4:]); q, r = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n), indexing="ij")
    X = q + a_ * np.sin(np.pi * q) * np.sin(2 * np.pi * r); Y = r + a_ * np.sin(2 * np.pi * q) * np.sin(np.pi * r)
elif f == "uniform":
    X, Y = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n), indexing="ij")
else:
    for a in csv.DictReader(open(f)): i, j = int(a["i"]) - 1, int(a["j"]) - 1; X[i, j] = float(a["x"]); Y[i, j] = float(a["y"])
    if Y[0, 0] > Y[0, -1]: X = X[:, ::-1]; Y = Y[:, ::-1]
    X[0, :] = 0; X[-1, :] = 1; Y[:, 0] = 0; Y[:, -1] = 1
o = [f"{n:5d}{n:5d}{NZ:5d}{0:5d}{0:5d}{0:5d}{0:5d}{0:5d}{0:5d}"]
for k in range(NZ):
    for i in range(n):
        o += ["", "", ""]
        for j in range(n):
            b = i in (0, n - 1) or j in (0, n - 1)
            t = X[i, j] if b else 0.5
            o.append(f"{X[i, j]:8.4f}{Y[i, j]:8.4f}{(k - 0.5) * dz:8.4f}{0:12.4E}{0:12.4E}{0:12.4E}  {0:11.4E}{t:11.4E}{0:11.4E}{0:11.4E}{0:11.4E}{0:11.4E}")
open("fort.8", "w").write("\n".join(o) + "\n")
open("in5", "w").write(" 3000    1   10    1\n 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000\n   20.0 0.0000\n")
