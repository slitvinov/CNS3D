import numpy as np, sys
L,M=(int(sys.argv[2]),int(sys.argv[3])) if len(sys.argv)>3 else (45,41)
a=np.array([[float(l[0:8]),float(l[8:16]),float(l[24:36])] for l in open(sys.argv[1]) if len(l)>100]).reshape(4,L,M,3)
x=a[1,:,0,0]; u=a[1,:,1,2]
i=next((i for i in range(2,L) if u[i-1]<0<=u[i]),None)
xr=np.interp(0,[u[i-1],u[i]],[x[i-1],x[i]]) if i else float("nan")
print(f"Xr/h = {xr:.2f}")
