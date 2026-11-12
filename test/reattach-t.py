import numpy as np, sys
L,M=45,(int(sys.argv[2]) if len(sys.argv)>2 else 43)
a=np.array([[float(l[0:8]),float(l[8:16]),float(l[24:36])] for l in open(sys.argv[1]) if len(l)>100]).reshape(4,L,M,3)
x=a[1,:,0,0]; u=a[1,:,1,2]
i=next((i for i in range(10,L) if u[i-1]<0<=u[i]),None)
print(f"Xr/h = {np.interp(0,[u[i-1],u[i]],[x[i-1],x[i]]):.2f}" if i else "no reattachment")
