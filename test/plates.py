import sys
L, M, NZ = 161, 21, 4
re = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
sk = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
x = [8.0 * i / (L - 1) for i in range(L)]; y = [j / (M - 1) for j in range(M)]; dz = y[1]
o = [f"{L:5d}{M:5d}{NZ:5d}{0:5d}{0:5d}{0:5d}{0:5d}{0:5d}{0:5d}"]
for k in range(NZ):
    for i in range(L):
        o += ["", "", ""]
        for j in range(M):
            wall = j == 0 or j == M - 1
            u = 0.0 if wall else 6 * y[j] * (1 - y[j])
            t = 0.0 if wall else 1.0
            o.append(f"{x[i] + sk * y[j]:8.4f}{y[j]:8.4f}{(k - 0.5) * dz:8.4f}{u:12.4E}{0:12.4E}{0:12.4E}  {0:11.4E}{t:11.4E}{0:11.4E}{0:11.4E}{0:11.4E}{0:11.4E}")
open("fort.8", "w").write("\n".join(o) + "\n")
open("in5", "w").write(" 3000    1   10    1\n 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000 0.5000\n" + f"{re:7.1f} 0.0000\n")
