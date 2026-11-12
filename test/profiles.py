import numpy as np, json, sys
L,M,N=20,17,9
a=np.array([[float(l[0:8]),float(l[8:16]),float(l[16:24]),float(l[24:36]),float(l[36:48]),float(l[48:60])] for l in open(sys.argv[1]) if len(l)>100]).reshape(N,L,M,6)
res={}
for nm,I in (("0",4),("60",10),("90",13)):
    t=np.radians(int(nm)); tx,ty=np.cos(t),np.sin(t)
    for zn,K in (("0",0),("0.25",4)):
        us=a[K,I,:,3]*tx+a[K,I,:,4]*ty
        rd=np.array([0.0, 0.0308, 0.0725, 0.1259, 0.1871, 0.2614, 0.3363, 0.4163, 0.4978, 0.5819, 0.6614, 0.7411, 0.8118, 0.8757, 0.9273, 0.9683, 1.0])[::-1]
        res[f"{nm}/{zn}"]=[list(rd),list(us)]
        print(f"theta={nm:>2} z/d={zn:4}: u/Uref from inner (r/d=0) to outer:", " ".join(f"{q:.2f}" for q in us[::-1]))
for I in (0,5,11,14,19):
    print("station",I+1,"flux ~", round(float(sum((a[K,I,:,3]**2+a[K,I,:,4]**2)**0.5 for K in range(N)).sum()),2))
json.dump(res,open("prof.json","w"))
