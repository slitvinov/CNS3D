import sys
n=int(sys.argv[1]); L=M=n; NZ=4; dz=1.0/(n-1)
o=[f"{L} {M} {NZ} 0 0 0 0 0 0"]
for k in range(NZ):
    for i in range(L):
        for j in range(M):
            u=1.0 if (j==M-1 and 0<i<L-1) else 0.0
            o.append(f"{i/(n-1):.6f} {j/(n-1):.6f} {(k-0.5)*dz:.6f} {u} 0 0 0 0")
open("field.txt","w").write("\n".join(o)+"\n")
