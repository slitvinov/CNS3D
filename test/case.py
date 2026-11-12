import sys
src, out, opts = sys.argv[1], sys.argv[2], sys.argv[3:]
s = open(src).read(); log = []
def rep(a, b, why, count=1):
    global s
    n = s.count(a)
    if n != count: sys.exit(f"edit not applicable ({n} matches): {a!r}")
    s = s.replace(a, b); log.append(why)
for o in opts:
    k, _, v = o.partition("=")
    if k == "dims":
        L0, M0, N0 = v.split(","); n = s.count("21,18,10")
        rep("21,18,10", f"{L0},{M0},{N0}", f"arrays ({L0},{M0},{N0})", n)
        mx = max(int(L0), int(M0), int(N0))
        rep("A(20),B(20),C(20),D(20)", f"A({mx}),B({mx}),C({mx}),D({mx})", "LINERX work arrays")
    elif k == "2d":
        rep("      IF(J.EQ.1.OR.J.EQ.M.OR.K.EQ.N) MC(I,J,K)=1\n", "      IF(J.EQ.1.OR.J.EQ.M) MC(I,J,K)=1\n", "no wall at K=N")
        rep("      DE(I,J,K)=DE(I,J,K+1)\n  10  CONTINUE\n",
            "      DE(I,J,K)=DE(I,J,K+1)\n      U(I,J,N)=U(I,J,N-1)\n      V(I,J,N)=V(I,J,N-1)\n      W(I,J,N)=0.0\n"
            "      TM(I,J,N)=TM(I,J,N-1)\n      DK(I,J,N)=DK(I,J,N-1)\n      DE(I,J,N)=DE(I,J,N-1)\n  10  CONTINUE\n", "K=N mirrors K=N-1")
    elif k == "walls":
        import re
        m = re.search(r"\n      IF\(([^)]*(?:\([^)]*\))*[^)]*)\) MC\(I,J,K\)=1\n      IF\(I\.GE\.L1", s)
        old = m.group(0); new = old.replace(m.group(1), v, 1)
        rep(old, new, f"walls where {v}")
    elif k == "outlet" and v == "u":
        rep("      UC=(V(I,J,K)+V(I,J-1,K)+V(I,J,K-1)+V(I,J-1,K-1))*0.25\n", "      UC=(U(I,J,K)+U(I,J-1,K)+U(I,J,K-1)+U(I,J-1,K-1))*0.25\n", "outlet flow from U")
        rep("      U(I+1,J,K)=U(I,J,K)\n      V(I+1,J,K)=V(I,J,K)-UINC\n", "      U(I+1,J,K)=U(I,J,K)-UINC\n      V(I+1,J,K)=V(I,J,K)\n", "outlet correction on U")
        if "      DO 50 J=2,JT+1\n" not in s:
            rep("      ARDEN=0.0\n      DO 50 J=2,JT\n      DO 50 K=2,KT\n", "      ARDEN=0.0\n      DO 50 J=2,JT+1\n      DO 50 K=2,KT+1\n", "outlet flow over the same cells as the inlet")
    elif k == "outletsum":
        rep("      ARDEN=0.0\n      DO 50 J=2,JT\n      DO 50 K=2,KT\n", "      ARDEN=0.0\n      DO 50 J=2,JT+1\n      DO 50 K=2,KT+1\n", "outlet flow over the same cells as the inlet")
    elif k == "tmsource":
        rep("  4   CONTINUE\n      GO TO 60\n",
            "  4   CONTINUE\n      DO 45 I=IS,IT\n      DO 45 J=JS,JT\n      DO 45 K=KS,KT\n"
            "      SU(I,J,K)=APO(I,J,K)*FO(I,J,K)*TJO(I,J,K)\n      SP(I,J,K)=-APO(I,J,K)*TJO(I,J,K)\n  45  CONTINUE\n      GO TO 60\n", "TM source: the time term only, without the mass imbalance of the cell")
        rep("      W(I+1,J,K)=W(I,J,K)\n      DK(I+1,J,K)=DK(I,J,K)\n      DE(I+1,J,K)=DE(I,J,K)\n  60  CONTINUE\n",
            "      W(I+1,J,K)=W(I,J,K)\n      TM(I+1,J,K)=TM(I,J,K)\n      DK(I+1,J,K)=DK(I,J,K)\n      DE(I+1,J,K)=DE(I,J,K)\n  60  CONTINUE\n", "TM carried to the outlet")
    elif k == "tmcross":
        rep("      SU(I,J,K)=APO(I,J,K)*FO(I,J,K)*TJO(I,J,K)\n      SP(I,J,K)=-APO(I,J,K)*TJO(I,J,K)\n  45  CONTINUE\n",
            "      GAE=0.5*(VISE(I+1,J,K)+VISE(I,J,K))/SIGF\n"
            "      GAW=0.5*(VISE(I-1,J,K)+VISE(I,J,K))/SIGF\n"
            "      GAN=0.5*(VISE(I,J+1,K)+VISE(I,J,K))/SIGF\n"
            "      GAS=0.5*(VISE(I,J-1,K)+VISE(I,J,K))/SIGF\n"
            "      GAT=0.5*(VISE(I,J,K+1)+VISE(I,J,K))/SIGF\n"
            "      GAB=0.5*(VISE(I,J,K-1)+VISE(I,J,K))/SIGF\n"
            "      DDNE=GAE*TYXE(I,J,K)+GAN*TXYN(I,J,K)\n"
            "      DDSE=-GAE*TYXE(I,J,K)-GAS*TXYS(I,J,K)\n"
            "      DDNW=-GAW*TYXW(I,J,K)-GAN*TXYN(I,J,K)\n"
            "      DDSW=GAW*TYXW(I,J,K)+GAS*TXYS(I,J,K)\n"
            "      DDET=GAE*TZXE(I,J,K)+GAT*TXZT(I,J,K)\n"
            "      DDEB=-GAE*TZXE(I,J,K)-GAB*TXZB(I,J,K)\n"
            "      DDWT=-GAW*TZXW(I,J,K)-GAT*TXZT(I,J,K)\n"
            "      DDWB=GAW*TZXW(I,J,K)+GAB*TXZB(I,J,K)\n"
            "      DDNT=GAN*TZYN(I,J,K)+GAT*TYZT(I,J,K)\n"
            "      DDNB=-GAN*TZYN(I,J,K)-GAB*TYZB(I,J,K)\n"
            "      DDST=-GAS*TZYS(I,J,K)-GAT*TYZT(I,J,K)\n"
            "      DDSB=GAS*TZYS(I,J,K)+GAB*TYZB(I,J,K)\n"
            "      SU(I,J,K)=DDNE*F(I+1,J+1,K)+DDSE*F(I+1,J-1,K)+\n"
            "     1         DDNW*F(I-1,J+1,K)+DDSW*F(I-1,J-1,K)+\n"
            "     2         DDET*F(I+1,J,K+1)+DDEB*F(I+1,J,K-1)+\n"
            "     3         DDWT*F(I-1,J,K+1)+DDWB*F(I-1,J,K-1)+\n"
            "     4         DDNT*F(I,J+1,K+1)+DDNB*F(I,J+1,K-1)+\n"
            "     5         DDST*F(I,J-1,K+1)+DDSB*F(I,J-1,K-1)\n"
            "      SU(I,J,K)=(SU(I,J,K)+APO(I,J,K)*FO(I,J,K))*TJO(I,J,K)\n"
            "      SP(I,J,K)=-APO(I,J,K)*TJO(I,J,K)\n  45  CONTINUE\n", "TM source: the cross-derivative diffusion of a skewed grid")
    elif k == "wallk":
        rep("      SU1=TJO1*(SUK1*DK1+(P1*TAUN1**2/VISE1+GENR))\n",
            "      SU1=SU1+TJO1*(P1*TAUN1**2/VISE1+GENR-GEN1)\n", "near-wall K: only the production replaced")
    elif k == "outlet" and v == "none":
        rep("      I=IT\n      FLOW=0.0\n", "      RETURN\n      I=IT\n      FLOW=0.0\n", "no outlet: I=L fixed")
        rep("      I=IT\n      DO 200 J=2,JT\n", "      RETURN\n      I=IT\n      DO 200 J=2,JT\n", "keep the east link")
    elif k == "bounc":
        for f, a in (("F(I,J-1,K)", "AS"), ("F(I+1,J,K)", "AE"), ("F(I-1,J,K)", "AW"), ("F(I,J,K+1)", "AT"), ("F(I,J,K-1)", "AB")):
            rep(f"SINZ(III),{f},AN(I,J,K)", f"SINZ(III),{f},{a}(I,J,K)", f"BOUNC wall link {a}")
    elif k == "unsteady":
        rep("     1     AB(I,J,K)+APO(I,J,K)\n      AP(I,J,K)=ANAB-SP(I,J,K)\n", "     1     AB(I,J,K)\n      AP(I,J,K)=ANAB-SP(I,J,K)\n", "unsteady term once")
    elif k == "heat":
        rep("      INSOT=.FALSE.\n", "      INSOT=.TRUE.\n", "temperature equation on")
    elif k == "wallvis":
        rep("      DO 22 I=IS,IT\n      DO 22 J=JS,JT\n      DO 22 K=KS,KT\n", "      DO 22 I=1,L\n      DO 22 J=1,M\n      DO 22 K=1,N\n", "F1 on all nodes")
    elif k == "erext":
        rep("      EREXT=1.E-3\n", f"      EREXT={v}\n", f"convergence tolerance {v}")
    elif k == "wallkt":
        rep("      P1=DDU**2+DDV**2+DDW**2\n", "      P1=(SINX1*DDU)**2+(SINY1*DDV)**2+(SINZ1*DDW)**2\n", "near-wall K production from the velocity along the wall")
    elif k == "pecap":
        for d in "EWNSTB":
            rep(f"AMAX1(ABS(0.5*C{d}),DD{d})", f"AMAX1(ABS(C{d}/{v}.),DD{d})", f"hybrid switch at cell Peclet {v} ({d})")
    elif k == "hyb":
        for d in "EWNSTB":
            rep(f"AMAX1(ABS(0.5*C{d}),DD{d})", f"AMAX1(ABS({v}*C{d}),DD{d})", f"hybrid factor {v} ({d})")
    elif k == "c1":
        rep("      C1=1.43\n", f"      C1={v}\n", f"C1 = {v}")
    elif k == "genshear":
        rep("      GEN(I,J,K)=VISE(I,J,K)*((UY+VX)**2+(VZ+WY)**2+(WX+UZ)**2+\n     1         2*(UX*UX+VY*VY+WZ*WZ))\n",
            "      GEN(I,J,K)=VISE(I,J,K)*((UY+VX)**2+(VZ+WY)**2+(WX+UZ)**2)\n", "k production from the shear terms only")
    elif k == "walldist":
        for a, b, d in (("      J1=J-1\n      J2=J-2\n", "      J1=J+1\n      J2=J-1\n", "N"), ("      J1=J+1\n      J2=J+2\n", "      J1=J-1\n      J2=J+1\n", "S"),
                        ("      I1=I-1\n      I2=I-2\n", "      I1=I+1\n      I2=I-1\n", "E"), ("      I1=I+1\n      I2=I+2\n", "      I1=I-1\n      I2=I+1\n", "W"),
                        ("      K1=K-1\n      K2=K-2\n", "      K1=K+1\n      K2=K-1\n", "T"), ("      K1=K+1\n      K2=K+2\n", "      K1=K-1\n      K2=K+1\n", "B")):
            rep(a, b, f"wall distance to the wall node ({d})")
    elif k == "nolip":
        rep("      IF(MC(I-1,J,K) .EQ. 0) GO TO 5\n", "      IF(MC(I-1,J,K) .EQ. 0) GO TO 5\n      IF(J.EQ.M2) GO TO 5\n", "no wall function at the node behind the step corner")
    elif k == "pe10":
        for c in "EWNSTB":
            rep(f"AMAX1(ABS(0.5*C{c}),DD{c})", f"AMAX1(ABS(0.1*C{c}),DD{c})", f"hybrid switch at cell Peclet 10 ({c})")
    else: sys.exit(f"unknown option {o}")
open(out, "w").write(s); print(out + ": " + "; ".join(log))
