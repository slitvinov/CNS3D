import sys, numpy as np
a = np.array([[float(l[0:8]), float(l[73:84])] for l in open(sys.argv[1]) if len(l) > 100]).reshape(4, 31, 31, 2)[1]
e = a[:, :, 1] - a[:, :, 0]; print(f"rms |T - x| {np.sqrt((e ** 2).mean()):.5f}, max {np.abs(e).max():.5f}")
