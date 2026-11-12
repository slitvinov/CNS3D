import numpy as np
L,M,N=20,17,9
up=[-2.7447, -1.8373, -1.031, -0.3901]; th=np.radians(np.arange(0,91,10)); dn=[0.1865, 0.4345, 0.739, 1.0906, 1.4892, 1.9318]
r=1.8+np.array([0.0, 0.0308, 0.0725, 0.1259, 0.1871, 0.2614, 0.3363, 0.4163, 0.4978, 0.5819, 0.6614, 0.7411, 0.8118, 0.8757, 0.9273, 0.9683, 1.0])[::-1]
z=np.array([0.0, 0.0588, 0.183, 0.3268, 0.4902, 0.6601, 0.8105, 0.9346, 1.0])*0.5
def prof(s,zz):
    Y=s-0.5; Z=zz; v=0*Y
    for n in range(1,201,2): v+=(-1)**((n-1)//2)*(1-np.cosh(n*np.pi*Z)/np.cosh(n*np.pi*0.5))*np.cos(n*np.pi*Y)/n**3
    return v
ss=np.linspace(0,1,401); zz=np.linspace(-0.5,0.5,401); S,Z=np.meshgrid(ss,zz); bulk=prof(S,Z).mean()
out=[f"{L} {M} {N} 0 0 0 0 0 0"]
for k in range(N):
    for i in range(L):
        for j in range(M):
            s=2.8-r[j]; q=prof(s,z[k])/bulk if 0<j<M-1 and k<N-1 else 0.0
            if i<4:   x,y,tx,ty=up[i],-r[j],1,0
            elif i<14: t=th[i-4]; x,y,tx,ty=r[j]*np.sin(t),-r[j]*np.cos(t),np.cos(t),np.sin(t)
            else:     x,y,tx,ty=r[j],dn[i-14],0,1
            out.append(f"{x:.6f} {y:.6f} {z[k]:.6f} {q*tx:.6e} {q*ty:.6e} 0.0 0.0 0.0")
open("field.txt","w").write("\n".join(out)+"\n")
