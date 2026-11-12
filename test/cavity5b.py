import csv, numpy as np
n = 31; NZ = 4; dz = 1.0 / (n - 1)
X = np.zeros((n, n)); Y = np.zeros((n, n))
for a in csv.DictReader(open("../redraw/fig5b_nodes.csv")):
    i, j = int(a["i"]) - 1, int(a["j"]) - 1
    X[i, j] = float(a["x"]); Y[i, j] = float(a["y"])
if Y[0, 0] > Y[0, -1]: X = X[:, ::-1]; Y = Y[:, ::-1]
X[0, :] = 0; X[-1, :] = 1; Y[:, 0] = 0; Y[:, -1] = 1
o = [f"{n} {n} {NZ} 0 0 0 0 0 0"]
for k in range(NZ):
    for i in range(n):
        for j in range(n):
            u = 1.0 if (j == n - 1 and 0 < i < n - 1) else 0.0
            o.append(f"{X[i, j]:.6f} {Y[i, j]:.6f} {(k - 0.5) * dz:.6f} {u} 0 0 0 0")
open("field.txt", "w").write("\n".join(o) + "\n")
