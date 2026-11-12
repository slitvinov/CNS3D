import sys, numpy as np
def ld(f): return np.array([[float(l[0:8]),float(l[8:16]),float(l[24:36]),float(l[36:48]),float(l[84:95])] for l in open(f) if len(l)>100])
a=ld(sys.argv[1]); b=ld(sys.argv[2])
for n,q in (("u",2),("v",3),("k",4)): print(n,"max |diff| %.2e (max %.2e)"%(np.abs(a[:,q]-b[:,q]).max(),np.abs(a[:,q]).max()))
for f,x in ((sys.argv[1],a),(sys.argv[2],b)):
    L=int(sys.argv[3]); A=x.reshape(4,L,int(sys.argv[4]),5)[1]; xx=A[:,0,0]; u=A[:,1,2]; i=next((i for i in range(10,L) if u[i-1]<0<=u[i]),None)
    print(f.split("/")[-1], "Xr=%.3f kmax=%.4f"%(np.interp(0,[u[i-1],u[i]],[xx[i-1],xx[i]]), A[xx>0,:,4].max()))
