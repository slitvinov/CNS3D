import numpy as np, sys
n=int(sys.argv[2]); a=np.array([[float(l[8:16]),float(l[24:36])] for l in open(sys.argv[1]) if len(l)>100]).reshape(4,n,n,2)
for y,u in a[1,n//2]: print(f"{y:.4f} {u:.4f}")
