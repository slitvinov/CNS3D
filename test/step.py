import numpy as np, sys
a=sys.argv[1:]
L,M,b=(int(a[0]),int(a[1]),float(a[2])) if len(a)>2 else (45,41,3.0)
ym,by,LX=(a[3],float(a[4]),float(a[5])) if len(a)>5 else ("uni",0.0,30.0); NZ=4
s=np.linspace(0,1,L); x=LX*(np.exp(b*s)-1)/(np.exp(b)-1) if b>0 else LX*s
t=np.linspace(-1,1,M)
if ym=="lip" and by>0: y=1+np.sign(t)*(1-np.tanh(by*(1-abs(t)))/np.tanh(by))
elif ym=="wall" and by>0: y=1+np.tanh(by*t)/np.tanh(by)
else: y=1+t
dz=2/(M-1)
o=[f"{L} {M} {NZ} 0 0 0 0 0 0"]
for k in range(NZ):
    for i in range(L):
        for j in range(M):
            e=y[j]-1; u=6*e*(1-e) if 0<e<1 else 0.0
            o.append(f"{x[i]:.6f} {y[j]:.6f} {(k-0.5)*dz:.6f} {u:.6e} 0 0 0 0")
open("field.txt","w").write("\n".join(o)+"\n"); print("dx at step %.3f, at exit %.3f"%(x[1]-x[0],x[-1]-x[-2]))
