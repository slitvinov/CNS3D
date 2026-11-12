      PROGRAM  CNS3D
C      BY: Y. S. CHEN       CN: 8/13/1985
C
C*****[ CURVILINEAR N-S CODE FOR 3-D INCOMPRESSIBLE FLOWS ]***********
C
      COMMON
     1/VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2  DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3  PP(21,18,10),W(21,18,10),TM(21,18,10)
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/PCCR/ DU(21,18,10),DV(21,18,10),DW(21,18,10)
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
     1/COEF/ AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2   SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3   AS(21,18,10),AT(21,18,10),AB(21,18,10),APO(21,18,10)
      COMMON
     1/TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJO(21,18,10),
     2   CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3   EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     3   SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
      COMMON
     1/TTRAN/TXXE(21,18,10),TXXW(21,18,10),TYYN(21,18,10),
     2   TYYS(21,18,10),TZZT(21,18,10),TZZB(21,18,10),
     3   TYXE(21,18,10),TYXW(21,18,10),TYZT(21,18,10),
     4   TYZB(21,18,10),TXYN(21,18,10),TXYS(21,18,10),
     5   TXZT(21,18,10),TXZB(21,18,10),TZXE(21,18,10),
     6   TZXW(21,18,10),TZYN(21,18,10),TZYS(21,18,10)
     1/UNSTDY/UO(21,18,10),VO(21,18,10),WO(21,18,10),DKO(21,18,10),
     2   DEO(21,18,10),DENO(21,18,10),TMO(21,18,10)
      LOGICAL INSOU,INSOV,INSOP,INSOK,INSOE,INPRO,INSOW,INSOT
C*****[ INPUT DATA GUIDE ]********************************************
C
C  NLIMT : MAXIMUM NO. OF ITERATIONS LIMIT
C
C  IG    = 1 : LAMINAR
C          2 : TURBULENT (K-E MODEL)
C
C  ISWP  : NO. OF SWEEPS FOR SOLVING THE P' EQUATION (PP).
C
C  ITT   : TOTAL NO. OF TIME STEPS.
C
C  ALU,ALV,ALW,ALP,ALK,ALE,ALVIS,ALC : UNDER-RELAXATION FACTORS
C
C  RENL  : REFERENCE REYNOLDS NUMBER.
C
C  DTT   : TIME STEP FOR UNSTEADY PROBLEMS.
C
C*********************************************************************
C-----INPUT DATA (PROBLEM CONTROL SETTING)
      READ(5,100) NLIMT,IG,ISWP,ITT
      READ(5,200) ALU,ALV,ALW,ALP,ALK,ALE,ALVIS,ALC
      READ(5,200) RENL,DTT
C-----CONSTANTS
      EREXT=1.E-3
      ISWU=7
      ISWV=7
      ISWW=7
      ISWK=5
      ISWE=5
      DENIN=1.0
      VISC=1./RENL
      SIGU=1.0
      SIGK=1.0
      SIGE=1.3
      CMU=0.09
      C1=1.43
      C2=1.92
      E=9.01069
      CK=0.4
      PI=3.141592654
      HINUM=1.E30
      SMNUM=1.E-30
      CMU1=CMU**0.25
      CMU2=CMU**0.75
      INSOU=.TRUE.
      INSOV=.TRUE.
      INSOP=.TRUE.
      INSOK=.TRUE.
      INSOE=.TRUE.
      INPRO=.TRUE.
      INSOW=.TRUE.
      INSOT=.FALSE.
      IF(IG .EQ. 2) GO TO 10
      INSOK=.FALSE.
      INSOE=.FALSE.
      INPRO=.FALSE.
  10  CONTINUE
C*****[READ IN INITIAL FLOW FIELDS FROM RESTART FILE (LU = 8)]********
      READ(8,100) L,M,N,L1,L2,M1,M2,N1,N2
      LO=L+1
      MO=M+1
      NO=N+1
      LT=L-1
      MT=M-1
      NT=N-1
C-----INITIALIZE VARIABLES
      CALL INIT
C-----RESTART FILE
      DO 50 K=1,N
      DO 50 I=1,L
      READ(8,400)
      DO 50 J=1,M
      READ(8,500) X(I,J,K),Y(I,J,K),Z(I,J,K),U(I,J,K),V(I,J,K),
     1 W(I,J,K),P(I,J,K),TM(I,J,K),DK(I,J,K),DE(I,J,K),
     2 VISETM,DENTM
      UO(I,J,K)=U(I,J,K)
      VO(I,J,K)=V(I,J,K)
      WO(I,J,K)=W(I,J,K)
      TMO(I,J,K)=TM(I,J,K)
      DKO(I,J,K)=DK(I,J,K)
      DEO(I,J,K)=DE(I,J,K)
      DENO(I,J,K)=DEN(I,J,K)
  50  CONTINUE
C*********************************************************************
C-----GET BOUNDARY CONTROL PARAMETERS
      CALL DIRCOS
C-----SET BOUNDARY TURBULENCE QUANTITIES TO ZERO
      DO 121 I=1,L
      DO 121 J=1,M
      DO 121 K=1,N
      IF(MC(I,J,K) .NE. 0) GO TO 122
      GO TO 121
 122  DK(I,J,K)=0.0
      DE(I,J,K)=0.0
      U(I,J,K)=0.0
      V(I,J,K)=0.0
      W(I,J,K)=0.0
 121  CONTINUE
C-----CALCULATE GRID TRANSFORMATION COEFFICIENTS
      CALL TRANF
C-----TURBULENT VISCOSITY
      IF(INPRO) CALL NEWVIS
C-----CALCULATE INLET MASS FLOW RATE
      FLOWIN=0.0
      I=1
      DO 45 J=2,M
      DO 45 K=2,N
      UC=(U(I,J,K)+U(I,J-1,K)+U(I,J,K-1)+U(I,J-1,K-1))*0.25
      DENC=(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I,J,K-1)+DEN(I,J-1,K-1))*0.25
      P1=(X(I,J,K)+X(I,J,K-1)-X(I,J-1,K)-X(I,J-1,K-1))*0.5
      P2=(Y(I,J,K)+Y(I,J,K-1)-Y(I,J-1,K)-Y(I,J-1,K-1))*0.5
      P3=(Z(I,J,K)+Z(I,J,K-1)-Z(I,J-1,K)-Z(I,J-1,K-1))*0.5
      Q1=(X(I,J,K)+X(I,J-1,K)-X(I,J,K-1)-X(I,J-1,K-1))*0.5
      Q2=(Y(I,J,K)+Y(I,J-1,K)-Y(I,J,K-1)-Y(I,J-1,K-1))*0.5
      Q3=(Z(I,J,K)+Z(I,J-1,K)-Z(I,J,K-1)-Z(I,J-1,K-1))*0.5
      AREA=SQRT(P1*P1+P2*P2+P3*P3)*SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      FLOWIN=FLOWIN+DENC*AREA*UC
  45  CONTINUE
      ITO=1
C-----TRANSIENT PROCESS
  2   CONTINUE
      CALL SYMOUT(3,1,2,L,2,M,2,N)
      ITER=1
C-----SOLUTION PROCEDURES START
  1   CONTINUE
      CALL SYMOUT(1,1,2,LT,2,MT,2,NT)
      IF(INSOU) CALL SOLVEQ(1,ISWU,ALU,SIGU,ERRU,U,UO)
      IF(INSOV) CALL SOLVEQ(2,ISWV,ALV,SIGU,ERRV,V,VO)
      IF(INSOW) CALL SOLVEQ(3,ISWW,ALW,SIGU,ERRW,W,WO)
      IF(INSOT) CALL SOLVEQ(4,ISWW,ALW,SIGU,ERRW,TM,TMO)
      IF(INSOK) CALL SOLVEQ(5,ISWK,ALK,SIGK,ERRK,DK,DKO)
      IF(INSOE) CALL SOLVEQ(6,ISWE,ALE,SIGE,ERRE,DE,DEO)
      IF(INSOP) CALL SOLVEQ(0,ISWP,ALP,SIGU,ERRM,PP,PP)
      IF(INPRO) CALL NEWVIS
C-----CONVERGENCE CHECK
      WRITE(6,300) ITER,ERRU,ERRV,ERRW,ERRM,ERRK,ERRE,U(7,2,6)
      ERRMAX=ERRM+ERRU+ERRV+ERRW
      IF(ITER .GE. 20 .AND. ERRMAX .GT. 1.E03) GO TO 99
      IF(ITER .GE. NLIMT .OR. ERRMAX .LE. EREXT) GO TO 99
      ITER=ITER+1
      GO TO 1
C-----PRINT OUT SOLUTIONS
  99  CONTINUE
      WRITE(7,100) L,M,N,L1,L2,M1,M2,N1,N2
      DO 901 K=1,N
      DO 901 I=1,L
      WRITE(7,400)
      DO 902 J=1,M
      XV=X(I,J,K)
      YV=Y(I,J,K)
      ZV=Z(I,J,K)
      WRITE(7,500) XV,YV,ZV,U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),
     1          TM(I,J,K),DK(I,J,K),DE(I,J,K),VISE(I,J,K),DEN(I,J,K)
  902  CONTINUE
  901  CONTINUE
      TIMT=DTT*ITO
      WRITE(7,300) ITO,TIMT
      IF(ITO .GE. ITT .OR. DTT .EQ. 0.0) GO TO 999
      ITO=ITO+1
      GO TO 2
 999  CONTINUE
 100  FORMAT(9I5)
 200  FORMAT(11F7.4)
 300  FORMAT(1X,I5,7E10.2)
 400  FORMAT(//)
 500  FORMAT(3F8.4,3E12.4,2X,6E11.4)
      STOP
      END
      SUBROUTINE DIRCOS
      COMMON
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
     1/TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJO(21,18,10),
     2   CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3   EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     3   SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
C-----SET DOMAIN BLOCKAGE CONTROL PARAMETER
C-----SCALAR BLOCKAGE :   MC(I,J,K)=1
C-----PRESSURE BLOCKAGE : MC(I,J,K)=2
      DO 10 I=1,LO
      DO 10 J=1,MO
      DO 10 K=1,NO
      MC(I,J,K)=0
      IF(J.EQ.1.OR.J.EQ.M.OR.K.EQ.N) MC(I,J,K)=1
      IF(I.GE.L1.AND.I.LE.L2.AND.J.GE.M1.AND.J.LE.M2.AND.
     1   K.GE.N1.AND.K.LE.N2) MC(I,J,K)=1
      IF(I.GT.L1.AND.I.LE.L2.AND.J.GT.M1.AND.J.LE.M2.AND.
     1   K.GT.N1.AND.K.LE.N2) MC(I,J,K)=2
C-----ADD BLOCKAGES AS NEEDED HERE
  10  CONTINUE
C-----CALCULATE BOUNDARY GRID SIZES AND ORIENTATIONS
      III=1
      DO 30 I=2,LT
      DO 30 J=2,MT
      DO 30 K=2,NT
      IF(MC(I,J,K) .NE. 0) GO TO 30
      MCT=MC(I+1,J,K)+MC(I-1,J,K)+MC(I,J+1,K)+MC(I,J-1,K)+
     1    MC(I,J,K+1)+MC(I,J,K-1)
      IF(MCT .EQ. 0) GO TO 30
      IF(MC(I,J+1,K) .EQ. 0) GO TO 2
C-----NORTH
      IBC(III)=I
      JBC(III)=J
      KBC(III)=K
      IITY(III)=1
      I1=I+1
      I2=I-1
      K1=K+1
      K2=K-1
      IF(I .EQ. L1) I2=I
      IF(I .EQ. L2) I1=I
      IF(K .EQ. N1) K2=K
      IF(K .EQ. N2) K1=K
      J1=J-1
      J2=J-2
      P1=(Y(I1,J,K1)-Y(I1,J,K2))*(Z(I1,J,K2)-Z(I2,J,K))-
     1   (Z(I1,J,K1)-Z(I1,J,K2))*(Y(I1,J,K2)-Y(I2,J,K))
      P2=(Z(I1,J,K1)-Z(I1,J,K2))*(X(I1,J,K2)-X(I2,J,K))-
     1   (X(I1,J,K1)-X(I1,J,K2))*(Z(I1,J,K2)-Z(I2,J,K))
      P3=(X(I1,J,K1)-X(I1,J,K2))*(Y(I1,J,K2)-Y(I2,J,K))-
     1   (Y(I1,J,K1)-Y(I1,J,K2))*(X(I1,J,K2)-X(I2,J,K))
      PQ=SQRT(P1*P1+P2*P2+P3*P3)
      P1=P1/PQ
      P2=P2/PQ
      P3=P3/PQ
      R1=(1.-P1**2)
      R2=(1.-P2**2)
      R3=(1.-P3**2)
      SINX(III)=SQRT(R1)
      SINY(III)=SQRT(R2)
      SINZ(III)=SQRT(R3)
      Q1=X(I,J,K)-X(I,J1,K)
      Q2=Y(I,J,K)-Y(I,J1,K)
      Q3=Z(I,J,K)-Z(I,J1,K)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      CC=1.0
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN(III)=BB*ABS(COTH)
      Q1=X(I,J,K)-X(I,J2,K)
      Q2=Y(I,J,K)-Y(I,J2,K)
      Q3=Z(I,J,K)-Z(I,J2,K)
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN1(III)=(BB*ABS(COTH)+YN(III))*0.5
      IJLO(I,J,K)=III
      III=III+1
  2   CONTINUE
      IF(MC(I,J-1,K) .EQ. 0) GO TO 3
C-----SOUTH
      IBC(III)=I
      JBC(III)=J
      KBC(III)=K
      IITY(III)=2
      I1=I+1
      I2=I-1
      K1=K+1
      K2=K-1
      J1=J+1
      J2=J+2
      P1=(Y(I1,J,K1)-Y(I1,J,K2))*(Z(I1,J,K2)-Z(I2,J,K))-
     1   (Z(I1,J,K1)-Z(I1,J,K2))*(Y(I1,J,K2)-Y(I2,J,K))
      P2=(Z(I1,J,K1)-Z(I1,J,K2))*(X(I1,J,K2)-X(I2,J,K))-
     1   (X(I1,J,K1)-X(I1,J,K2))*(Z(I1,J,K2)-Z(I2,J,K))
      P3=(X(I1,J,K1)-X(I1,J,K2))*(Y(I1,J,K2)-Y(I2,J,K))-
     1   (Y(I1,J,K1)-Y(I1,J,K2))*(X(I1,J,K2)-X(I2,J,K))
      PQ=SQRT(P1*P1+P2*P2+P3*P3)
      P1=P1/PQ
      P2=P2/PQ
      P3=P3/PQ
      R1=(1.-P1**2)
      R2=(1.-P2**2)
      R3=(1.-P3**2)
      SINX(III)=SQRT(R1)
      SINY(III)=SQRT(R2)
      SINZ(III)=SQRT(R3)
      Q1=X(I,J,K)-X(I,J1,K)
      Q2=Y(I,J,K)-Y(I,J1,K)
      Q3=Z(I,J,K)-Z(I,J1,K)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      CC=1.0
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN(III)=BB*ABS(COTH)
      Q1=X(I,J,K)-X(I,J2,K)
      Q2=Y(I,J,K)-Y(I,J2,K)
      Q3=Z(I,J,K)-Z(I,J2,K)
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN1(III)=(BB*ABS(COTH)+YN(III))*0.5
      IJLO(I,J,K)=III
      III=III+1
  3   CONTINUE
      IF(MC(I+1,J,K) .EQ. 0) GO TO 4
C-----EAST
      IBC(III)=I
      JBC(III)=J
      KBC(III)=K
      IITY(III)=3
      J1=J+1
      J2=J-1
      K1=K+1
      K2=K-1
      IF(J .EQ. M1) J2=J
      IF(J .EQ. M2) J1=J
      IF(K .EQ. N1) K2=K
      IF(K .EQ. N2) K1=K
      I1=I-1
      I2=I-2
      P1=(Y(I,J1,K1)-Y(I,J1,K2))*(Z(I,J1,K2)-Z(I,J2,K))-
     1   (Z(I,J1,K1)-Z(I,J1,K2))*(Y(I,J1,K2)-Y(I,J2,K))
      P2=(Z(I,J1,K1)-Z(I,J1,K2))*(X(I,J1,K2)-X(I,J2,K))-
     1   (X(I,J1,K1)-X(I,J1,K2))*(Z(I,J1,K2)-Z(I,J2,K))
      P3=(X(I,J1,K1)-X(I,J1,K2))*(Y(I,J1,K2)-Y(I,J2,K))-
     1   (Y(I,J1,K1)-Y(I,J1,K2))*(X(I,J1,K2)-X(I,J2,K))
      PQ=SQRT(P1*P1+P2*P2+P3*P3)
      P1=P1/PQ
      P2=P2/PQ
      P3=P3/PQ
      R1=(1.-P1**2)
      R2=(1.-P2**2)
      R3=(1.-P3**2)
      SINX(III)=SQRT(R1)
      SINY(III)=SQRT(R2)
      SINZ(III)=SQRT(R3)
      Q1=X(I,J,K)-X(I1,J,K)
      Q2=Y(I,J,K)-Y(I1,J,K)
      Q3=Z(I,J,K)-Z(I1,J,K)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      CC=1.0
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN(III)=BB*ABS(COTH)
      Q1=X(I,J,K)-X(I2,J,K)
      Q2=Y(I,J,K)-Y(I2,J,K)
      Q3=Z(I,J,K)-Z(I2,J,K)
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN1(III)=(BB*ABS(COTH)+YN(III))*0.5
      IJLO(I,J,K)=III
      III=III+1
  4   CONTINUE
      IF(MC(I-1,J,K) .EQ. 0) GO TO 5
C-----WEST
      IBC(III)=I
      JBC(III)=J
      KBC(III)=K
      IITY(III)=4
      J1=J+1
      J2=J-1
      K1=K+1
      K2=K-1
      IF(J .EQ. M1) J2=J
      IF(J .EQ. M2) J1=J
      IF(K .EQ. N1) K2=K
      IF(K .EQ. N2) K1=K
      I1=I+1
      I2=I+2
      P1=(Y(I,J1,K1)-Y(I,J1,K2))*(Z(I,J1,K2)-Z(I,J2,K))-
     1   (Z(I,J1,K1)-Z(I,J1,K2))*(Y(I,J1,K2)-Y(I,J2,K))
      P2=(Z(I,J1,K1)-Z(I,J1,K2))*(X(I,J1,K2)-X(I,J2,K))-
     1   (X(I,J1,K1)-X(I,J1,K2))*(Z(I,J1,K2)-Z(I,J2,K))
      P3=(X(I,J1,K1)-X(I,J1,K2))*(Y(I,J1,K2)-Y(I,J2,K))-
     1   (Y(I,J1,K1)-Y(I,J1,K2))*(X(I,J1,K2)-X(I,J2,K))
      PQ=SQRT(P1*P1+P2*P2+P3*P3)
      P1=P1/PQ
      P2=P2/PQ
      P3=P3/PQ
      R1=1.-P1**2
      R2=1.-P2**2
      R3=1.-P3**2
      SINX(III)=SQRT(R1)
      SINY(III)=SQRT(R2)
      SINZ(III)=SQRT(R3)
      Q1=X(I,J,K)-X(I1,J,K)
      Q2=Y(I,J,K)-Y(I1,J,K)
      Q3=Z(I,J,K)-Z(I1,J,K)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      CC=1.0
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN(III)=BB*ABS(COTH)
      Q1=X(I,J,K)-X(I2,J,K)
      Q2=Y(I,J,K)-Y(I2,J,K)
      Q3=Z(I,J,K)-Z(I2,J,K)
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN1(III)=(BB*ABS(COTH)+YN(III))*0.5
      IJLO(I,J,K)=III
      III=III+1
  5   CONTINUE
      IF(MC(I,J,K+1) .EQ. 0) GO TO 6
C-----TOP
      IBC(III)=I
      JBC(III)=J
      KBC(III)=K
      IITY(III)=5
      I1=I+1
      I2=I-1
      J1=J+1
      J2=J-1
      IF(I .EQ. L1) I2=I
      IF(I .EQ. L2) I1=I
      IF(J .EQ. M1) J2=J
      IF(J .EQ. M2) J1=J
      K1=K-1
      K2=K-2
      P1=(Y(I1,J1,K)-Y(I1,J2,K))*(Z(I1,J2,K)-Z(I2,J,K))-
     1   (Z(I1,J1,K)-Z(I1,J2,K))*(Y(I1,J2,K)-Y(I2,J,K))
      P2=(Z(I1,J1,K)-Z(I1,J2,K))*(X(I1,J2,K)-X(I2,J,K))-
     1   (X(I1,J1,K)-X(I1,J2,K))*(Z(I1,J2,K)-Z(I2,J,K))
      P3=(X(I1,J1,K)-X(I1,J2,K))*(Y(I1,J2,K)-Y(I2,J,K))-
     1   (Y(I1,J1,K)-Y(I1,J2,K))*(X(I1,J2,K)-X(I2,J,K))
      PQ=SQRT(P1*P1+P2*P2+P3*P3)
      P1=P1/PQ
      P2=P2/PQ
      P3=P3/PQ
      R1=1.-P1**2
      R2=1.-P2**2
      R3=1.-P3**2
      SINX(III)=SQRT(R1)
      SINY(III)=SQRT(R2)
      SINZ(III)=SQRT(R3)
      Q1=X(I,J,K)-X(I,J,K1)
      Q2=Y(I,J,K)-Y(I,J,K1)
      Q3=Z(I,J,K)-Z(I,J,K1)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      CC=1.0
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN(III)=BB*ABS(COTH)
      Q1=X(I,J,K)-X(I,J,K2)
      Q2=Y(I,J,K)-Y(I,J,K2)
      Q3=Z(I,J,K)-Z(I,J,K2)
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN1(III)=(BB*ABS(COTH)+YN(III))*0.5
      IJLO(I,J,K)=III
      III=III+1
  6   CONTINUE
      IF(MC(I,J,K-1) .EQ. 0) GO TO 30
C-----BOTTOM
      IBC(III)=I
      JBC(III)=J
      KBC(III)=K
      IITY(III)=6
      I1=I+1
      I2=I-1
      J1=J+1
      J2=J-1
      IF(I .EQ. L1) I2=I
      IF(I .EQ. L2) I1=I
      IF(J .EQ. M1) J2=J
      IF(J .EQ. M2) J1=J
      K1=K+1
      K2=K+2
      P1=(Y(I1,J1,K)-Y(I1,J2,K))*(Z(I1,J2,K)-Z(I2,J,K))-
     1   (Z(I1,J1,K)-Z(I1,J2,K))*(Y(I1,J2,K)-Y(I2,J,K))
      P2=(Z(I1,J1,K)-Z(I1,J2,K))*(X(I1,J2,K)-X(I2,J,K))-
     1   (X(I1,J1,K)-X(I1,J2,K))*(Z(I1,J2,K)-Z(I2,J,K))
      P3=(X(I1,J1,K)-X(I1,J2,K))*(Y(I1,J2,K)-Y(I2,J,K))-
     1   (Y(I1,J1,K)-Y(I1,J2,K))*(X(I1,J2,K)-X(I2,J,K))
      PQ=SQRT(P1*P1+P2*P2+P3*P3)
      P1=P1/PQ
      P2=P2/PQ
      P3=P3/PQ
      R1=1.-P1**2
      R2=1.-P2**2
      R3=1.-P3**2
      SINX(III)=SQRT(R1)
      SINY(III)=SQRT(R2)
      SINZ(III)=SQRT(R3)
      Q1=X(I,J,K)-X(I,J,K1)
      Q2=Y(I,J,K)-Y(I,J,K1)
      Q3=Z(I,J,K)-Z(I,J,K1)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      CC=1.0
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN(III)=BB*ABS(COTH)
      Q1=X(I,J,K)-X(I,J,K2)
      Q2=Y(I,J,K)-Y(I,J,K2)
      Q3=Z(I,J,K)-Z(I,J,K2)
      BB=SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      AA=SQRT((Q1-P1)**2+(Q2-P2)**2+(Q3-P3)**2)
      COTH=(BB*BB+CC*CC-AA*AA)/(2*BB*CC)
      YN1(III)=(BB*ABS(COTH)+YN(III))*0.5
      IJLO(I,J,K)=III
      III=III+1
  30  CONTINUE
      IITO=III-1
      WRITE(6,100) LO,MO,NO,IITO
 100  FORMAT(4I5)
      RETURN
      END
      SUBROUTINE TRANF
      COMMON
     1/VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2  DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3  PP(21,18,10),W(21,18,10),TM(21,18,10)
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/PCCR/ DU(21,18,10),DV(21,18,10),DW(21,18,10)
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
     1/COEF/ AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2   SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3   AS(21,18,10),AT(21,18,10),AB(21,18,10),APO(21,18,10)
      COMMON
     1/TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJO(21,18,10),
     2   CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3   EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     3   SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
      COMMON
     1/TTRAN/TXXE(21,18,10),TXXW(21,18,10),TYYN(21,18,10),
     2   TYYS(21,18,10),TZZT(21,18,10),TZZB(21,18,10),
     3   TYXE(21,18,10),TYXW(21,18,10),TYZT(21,18,10),
     4   TYZB(21,18,10),TXYN(21,18,10),TXYS(21,18,10),
     5   TXZT(21,18,10),TXZB(21,18,10),TZXE(21,18,10),
     6   TZXW(21,18,10),TZYN(21,18,10),TZYS(21,18,10)
C-----CALCULATE GRID TRANSFORMATION COEFFICIENTS
      DO 40 I=2,L
      DO 40 J=2,M
      DO 40 K=2,N
      P1=(X(I,J,K)+X(I,J-1,K)+X(I,J,K-1)+X(I,J-1,K-1)-X(I-1,J,K)-
     1  X(I-1,J-1,K)-X(I-1,J,K-1)-X(I-1,J-1,K-1))*0.25
      P2=(X(I,J,K)+X(I-1,J,K)+X(I,J,K-1)+X(I-1,J,K-1)-X(I,J-1,K)-
     1  X(I-1,J-1,K)-X(I,J-1,K-1)-X(I-1,J-1,K-1))*0.25
      P3=(X(I,J,K)+X(I-1,J,K)+X(I,J-1,K)+X(I-1,J-1,K)-X(I,J,K-1)-
     1  X(I-1,J,K-1)-X(I,J-1,K-1)-X(I-1,J-1,K-1))*0.25
      Q1=(Y(I,J,K)+Y(I,J-1,K)+Y(I,J,K-1)+Y(I,J-1,K-1)-Y(I-1,J,K)-
     1  Y(I-1,J-1,K)-Y(I-1,J,K-1)-Y(I-1,J-1,K-1))*0.25
      Q2=(Y(I,J,K)+Y(I-1,J,K)+Y(I,J,K-1)+Y(I-1,J,K-1)-Y(I,J-1,K)-
     1  Y(I-1,J-1,K)-Y(I,J-1,K-1)-Y(I-1,J-1,K-1))*0.25
      Q3=(Y(I,J,K)+Y(I-1,J,K)+Y(I,J-1,K)+Y(I-1,J-1,K)-Y(I,J,K-1)-
     1  Y(I-1,J,K-1)-Y(I,J-1,K-1)-Y(I-1,J-1,K-1))*0.25
      R1=(Z(I,J,K)+Z(I,J-1,K)+Z(I,J,K-1)+Z(I,J-1,K-1)-Z(I-1,J,K)-
     1  Z(I-1,J-1,K)-Z(I-1,J,K-1)-Z(I-1,J-1,K-1))*0.25
      R2=(Z(I,J,K)+Z(I-1,J,K)+Z(I,J,K-1)+Z(I-1,J,K-1)-Z(I,J-1,K)-
     1  Z(I-1,J-1,K)-Z(I,J-1,K-1)-Z(I-1,J-1,K-1))*0.25
      R3=(Z(I,J,K)+Z(I-1,J,K)+Z(I,J-1,K)+Z(I-1,J-1,K)-Z(I,J,K-1)-
     1  Z(I-1,J,K-1)-Z(I,J-1,K-1)-Z(I-1,J-1,K-1))*0.25
      PTR=1.0/(P1*(Q2*R3-Q3*R2)-P2*(Q1*R3-Q3*R1)+P3*(Q1*R2-Q2*R1))
      SU(I,J,K)=ABS(1.0/PTR)
      AE(I,J,K)=PTR*(Q2*R3-Q3*R2)
      AW(I,J,K)=-PTR*(P2*R3-P3*R2)
      AN(I,J,K)=PTR*(P2*Q3-P3*Q2)
      AS(I,J,K)=-PTR*(Q1*R3-Q3*R1)
      AT(I,J,K)=PTR*(P1*R3-P3*R1)
      AB(I,J,K)=-PTR*(P1*Q3-P3*Q1)
      DU(I,J,K)=PTR*(Q1*R2-Q2*R1)
      DV(I,J,K)=-PTR*(P1*R2-P2*R1)
      DW(I,J,K)=PTR*(P1*Q2-P2*Q1)
  40  CONTINUE
      CALL WALVAL(1.0,2,L,2,M,2,N,AE)
      CALL WALVAL(1.0,2,L,2,M,2,N,AW)
      CALL WALVAL(1.0,2,L,2,M,2,N,AN)
      CALL WALVAL(1.0,2,L,2,M,2,N,AS)
      CALL WALVAL(1.0,2,L,2,M,2,N,AT)
      CALL WALVAL(1.0,2,L,2,M,2,N,AB)
      CALL WALVAL(1.0,2,L,2,M,2,N,DU)
      CALL WALVAL(1.0,2,L,2,M,2,N,DV)
      CALL WALVAL(1.0,2,L,2,M,2,N,DW)
      CALL WALVAL(1.0,2,L,2,M,2,N,SU)
      DO 80 I=1,L
      DO 80 J=1,M
      DO 80 K=1,N
      CX(I,J,K)=(AE(I,J,K)+AE(I,J+1,K)+AE(I,J,K+1)+AE(I,J+1,K+1)+
     1 AE(I+1,J,K)+AE(I+1,J+1,K)+AE(I+1,J,K+1)+AE(I+1,J+1,K+1))*0.125
      CY(I,J,K)=(AW(I,J,K)+AW(I,J+1,K)+AW(I,J,K+1)+AW(I,J+1,K+1)+
     1 AW(I+1,J,K)+AW(I+1,J+1,K)+AW(I+1,J,K+1)+AW(I+1,J+1,K+1))*0.125
      CZ(I,J,K)=(AN(I,J,K)+AN(I,J+1,K)+AN(I,J,K+1)+AN(I,J+1,K+1)+
     1 AN(I+1,J,K)+AN(I+1,J+1,K)+AN(I+1,J,K+1)+AN(I+1,J+1,K+1))*0.125
      EX(I,J,K)=(AS(I,J,K)+AS(I,J+1,K)+AS(I,J,K+1)+AS(I,J+1,K+1)+
     1 AS(I+1,J,K)+AS(I+1,J+1,K)+AS(I+1,J,K+1)+AS(I+1,J+1,K+1))*0.125
      EY(I,J,K)=(AT(I,J,K)+AT(I,J+1,K)+AT(I,J,K+1)+AT(I,J+1,K+1)+
     1 AT(I+1,J,K)+AT(I+1,J+1,K)+AT(I+1,J,K+1)+AT(I+1,J+1,K+1))*0.125
      EZ(I,J,K)=(AB(I,J,K)+AB(I,J+1,K)+AB(I,J,K+1)+AB(I,J+1,K+1)+
     1 AB(I+1,J,K)+AB(I+1,J+1,K)+AB(I+1,J,K+1)+AB(I+1,J+1,K+1))*0.125
      SX(I,J,K)=(DU(I,J,K)+DU(I,J+1,K)+DU(I,J,K+1)+DU(I,J+1,K+1)+
     1 DU(I+1,J,K)+DU(I+1,J+1,K)+DU(I+1,J,K+1)+DU(I+1,J+1,K+1))*0.125
      SY(I,J,K)=(DV(I,J,K)+DV(I,J+1,K)+DV(I,J,K+1)+DV(I,J+1,K+1)+
     1 DV(I+1,J,K)+DV(I+1,J+1,K)+DV(I+1,J,K+1)+DV(I+1,J+1,K+1))*0.125
      SZ(I,J,K)=(DW(I,J,K)+DW(I,J+1,K)+DW(I,J,K+1)+DW(I,J+1,K+1)+
     1 DW(I+1,J,K)+DW(I+1,J+1,K)+DW(I+1,J,K+1)+DW(I+1,J+1,K+1))*0.125
      TJO(I,J,K)=(SU(I,J,K)+SU(I,J+1,K)+SU(I,J,K+1)+SU(I,J+1,K+1)+
     1 SU(I+1,J,K)+SU(I+1,J+1,K)+SU(I+1,J,K+1)+SU(I+1,J+1,K+1))*0.125
  80  CONTINUE
      DO 200 I=1,LO
      DO 200 J=1,MO
      DO 200 K=1,NO
      DU(I,J,K)=0.0
 200  CONTINUE
      DO 160 I=2,LT
      DO 160 J=2,MT
      DO 160 K=2,NT
      CXE=(CX(I+1,J,K)+CX(I,J,K))*0.5
      CXW=(CX(I-1,J,K)+CX(I,J,K))*0.5
      CXN=(CX(I,J+1,K)+CX(I,J,K))*0.5
      CXS=(CX(I,J,K)+CX(I,J-1,K))*0.5
      CXT=(CX(I,J,K+1)+CX(I,J,K))*0.5
      CXB=(CX(I,J,K)+CX(I,J,K-1))*0.5
      CYE=(CY(I+1,J,K)+CY(I,J,K))*0.5
      CYW=(CY(I,J,K)+CY(I-1,J,K))*0.5
      CYN=(CY(I,J+1,K)+CY(I,J,K))*0.5
      CYS=(CY(I,J,K)+CY(I,J-1,K))*0.5
      CYT=(CY(I,J,K+1)+CY(I,J,K))*0.5
      CYB=(CY(I,J,K)+CY(I,J,K-1))*0.5
      CZE=(CZ(I+1,J,K)+CZ(I,J,K))*0.5
      CZW=(CZ(I,J,K)+CZ(I-1,J,K))*0.5
      CZN=(CZ(I,J+1,K)+CZ(I,J,K))*0.5
      CZS=(CZ(I,J,K)+CZ(I,J-1,K))*0.5
      CZT=(CZ(I,J,K+1)+CZ(I,J,K))*0.5
      CZB=(CZ(I,J,K)+CZ(I,J,K-1))*0.5
      EXE=(EX(I+1,J,K)+EX(I,J,K))*0.5
      EXW=(EX(I,J,K)+EX(I-1,J,K))*0.5
      EXN=(EX(I,J+1,K)+EX(I,J,K))*0.5
      EXS=(EX(I,J,K)+EX(I,J-1,K))*0.5
      EXT=(EX(I,J,K+1)+EX(I,J,K))*0.5
      EXB=(EX(I,J,K)+EX(I,J,K-1))*0.5
      EYE=(EY(I+1,J,K)+EY(I,J,K))*0.5
      EYW=(EY(I,J,K)+EY(I-1,J,K))*0.5
      EYN=(EY(I,J+1,K)+EY(I,J,K))*0.5
      EYS=(EY(I,J,K)+EY(I,J-1,K))*0.5
      EYT=(EY(I,J,K+1)+EY(I,J,K))*0.5
      EYB=(EY(I,J,K)+EY(I,J,K-1))*0.5
      EZE=(EZ(I+1,J,K)+EZ(I,J,K))*0.5
      EZW=(EZ(I,J,K)+EZ(I-1,J,K))*0.5
      EZN=(EZ(I,J+1,K)+EZ(I,J,K))*0.5
      EZS=(EZ(I,J,K)+EZ(I,J-1,K))*0.5
      EZT=(EZ(I,J,K+1)+EZ(I,J,K))*0.5
      EZB=(EZ(I,J,K)+EZ(I,J,K-1))*0.5
      SXE=(SX(I+1,J,K)+SX(I,J,K))*0.5
      SXW=(SX(I,J,K)+SX(I-1,J,K))*0.5
      SXN=(SX(I,J+1,K)+SX(I,J,K))*0.5
      SXS=(SX(I,J,K)+SX(I,J-1,K))*0.5
      SXT=(SX(I,J,K+1)+SX(I,J,K))*0.5
      SXB=(SX(I,J,K)+SX(I,J,K-1))*0.5
      SYE=(SY(I+1,J,K)+SY(I,J,K))*0.5
      SYW=(SY(I,J,K)+SY(I-1,J,K))*0.5
      SYN=(SY(I,J+1,K)+SY(I,J,K))*0.5
      SYS=(SY(I,J,K)+SY(I,J-1,K))*0.5
      SYT=(SY(I,J,K+1)+SY(I,J,K))*0.5
      SYB=(SY(I,J,K)+SY(I,J,K-1))*0.5
      SZE=(SZ(I+1,J,K)+SZ(I,J,K))*0.5
      SZW=(SZ(I,J,K)+SZ(I-1,J,K))*0.5
      SZN=(SZ(I,J+1,K)+SZ(I,J,K))*0.5
      SZS=(SZ(I,J,K)+SZ(I,J-1,K))*0.5
      SZT=(SZ(I,J,K+1)+SZ(I,J,K))*0.5
      SZB=(SZ(I,J,K)+SZ(I,J,K-1))*0.5
      TXXE(I,J,K)=CXE*CX(I,J,K)+CYE*CY(I,J,K)+CZE*CZ(I,J,K)
      TXXW(I,J,K)=CXW*CX(I,J,K)+CYW*CY(I,J,K)+CZW*CZ(I,J,K)
      TYYN(I,J,K)=EXN*EX(I,J,K)+EYN*EY(I,J,K)+EZN*EZ(I,J,K)
      TYYS(I,J,K)=EXS*EX(I,J,K)+EYS*EY(I,J,K)+EZS*EZ(I,J,K)
      TZZT(I,J,K)=SXT*SX(I,J,K)+SYT*SY(I,J,K)+SZT*SZ(I,J,K)
      TZZB(I,J,K)=SXB*SX(I,J,K)+SYB*SY(I,J,K)+SZB*SZ(I,J,K)
      TYXE(I,J,K)=(EXE*CX(I,J,K)+EYE*CY(I,J,K)+EZE*CZ(I,J,K))*0.25
      TYXW(I,J,K)=(EXW*CX(I,J,K)+EYW*CY(I,J,K)+EZW*CZ(I,J,K))*0.25
      TYZT(I,J,K)=(EXT*SX(I,J,K)+EYT*SY(I,J,K)+EZT*SZ(I,J,K))*0.25
      TYZB(I,J,K)=(EXB*SX(I,J,K)+EYB*SY(I,J,K)+EZB*SZ(I,J,K))*0.25
      TXYN(I,J,K)=(CXN*EX(I,J,K)+CYN*EY(I,J,K)+CZN*EZ(I,J,K))*0.25
      TXYS(I,J,K)=(CXS*EX(I,J,K)+CYS*EY(I,J,K)+CZS*EZ(I,J,K))*0.25
      TXZT(I,J,K)=(CXT*SX(I,J,K)+CYT*SY(I,J,K)+CZT*SZ(I,J,K))*0.25
      TXZB(I,J,K)=(CXB*SX(I,J,K)+CYB*SY(I,J,K)+CZB*SZ(I,J,K))*0.25
      TZXE(I,J,K)=(SXE*CX(I,J,K)+SYE*CY(I,J,K)+SZE*CZ(I,J,K))*0.25
      TZXW(I,J,K)=(SXW*CX(I,J,K)+SYW*CY(I,J,K)+SZW*CZ(I,J,K))*0.25
      TZYN(I,J,K)=(SXN*EX(I,J,K)+SYN*EY(I,J,K)+SZN*EZ(I,J,K))*0.25
      TZYS(I,J,K)=(SXS*EX(I,J,K)+SYS*EY(I,J,K)+SZS*EZ(I,J,K))*0.25
 160  CONTINUE
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TXXE)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TXXW)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TYYN)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TYYS)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TZZT)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TZZB)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TYXE)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TYXW)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TYZT)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TYZB)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TXYN)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TXYS)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TXZT)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TXZB)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TZXE)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TZXW)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TZYN)
      CALL WALVAL(1.0,2,LT,2,MT,2,NT,TZYS)
      RETURN
      END
      SUBROUTINE INIT
      COMMON
     1/VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2  DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3  PP(21,18,10),W(21,18,10),TM(21,18,10)
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/PCCR/ DU(21,18,10),DV(21,18,10),DW(21,18,10)
      COMMON
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
C-----INITIALIZE VARIABLES
      DO 10 I=1,LO
      DO 10 J=1,MO
      DO 10 K=1,NO
      U(I,J,K)=0.0
      V(I,J,K)=0.0
      W(I,J,K)=0.0
      P(I,J,K)=0.0
      PP(I,J,K)=0.0
      DK(I,J,K)=0.0
      DE(I,J,K)=0.0
      DEN(I,J,K)=DENIN
      VISE(I,J,K)=VISC
  10  CONTINUE
      RETURN
      END
      SUBROUTINE NEWVIS
      COMMON
     1/VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2  DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3  PP(21,18,10),W(21,18,10),TM(21,18,10)
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/PCCR/ DU(21,18,10),DV(21,18,10),DW(21,18,10)
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
      COMMON
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
C-----EVALUATE TURBULENT VISCOSITY
      DO 10 I=1,L
      DO 10 J=1,M
      DO 10 K=1,N
      IF(DK(I,J,K) .LE. SMNUM) DK(I,J,K)=SMNUM
      IF(DE(I,J,K) .LE. SMNUM) DE(I,J,K)=SMNUM
      IF(DE(I,J,K) .LE. SMNUM) GO TO 12
      TURVIS=DEN(I,J,K)*CMU*DK(I,J,K)**2/DE(I,J,K)+VISC
      GO TO 14
  12  TURVIS=VISC
  14  CONTINUE
      VISE(I,J,K)=VISE(I,J,K)+ALVIS*(TURVIS-VISE(I,J,K))
  10  CONTINUE
      RETURN
      END
      SUBROUTINE SOLVEQ(IE,ISWF,ALF,SIGF,ERRF,F,FO)
      DIMENSION F(21,18,10),F1(21,18,10),FO(21,18,10)
      COMMON
     1/VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2  DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3  PP(21,18,10),W(21,18,10),TM(21,18,10)
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/PCCR/ DU(21,18,10),DV(21,18,10),DW(21,18,10)
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
     1/COEF/ AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2   SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3   AS(21,18,10),AT(21,18,10),AB(21,18,10),APO(21,18,10)
      COMMON
     1/TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJO(21,18,10),
     2   CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3   EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     3   SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
      COMMON
     1/TTRAN/TXXE(21,18,10),TXXW(21,18,10),TYYN(21,18,10),
     2   TYYS(21,18,10),TZZT(21,18,10),TZZB(21,18,10),
     3   TYXE(21,18,10),TYXW(21,18,10),TYZT(21,18,10),
     4   TYZB(21,18,10),TXYN(21,18,10),TXYS(21,18,10),
     5   TXZT(21,18,10),TXZB(21,18,10),TZXE(21,18,10),
     6   TZXW(21,18,10),TZYN(21,18,10),TZYS(21,18,10)
C------TRANSPORT EQUATIONS LINERAIZATION AND SOLVER
      ERRF=0.0
      PI=3.141592654
C-----PRESSURE CORRECTION SOLVER STARTS FROM 10
      IF(IE .EQ. 0) GO TO 10
C-----U, V, W, TM, K & E EQUATIONS
      IS=2
      IT=LT
      JS=2
      JT=MT
      KS=2
      KT=NT
      GO TO (21,29,29,29,29,21,21,21), IE
  21  CONTINUE
      DO 22 I=1,L
      DO 22 J=1,M
      DO 22 K=1,N
  22  F1(I,J,K)=VISE(I,J,K)/SIGF
C-----EVALUATE LINK COEFF. AND SOURCE TERMS
      DO 20 I=IS,IT
      DO 20 J=JS,JT
      DO 20 K=KS,KT
      GAE=0.5*(F1(I+1,J,K)+F1(I,J,K))
      GAW=0.5*(F1(I-1,J,K)+F1(I,J,K))
      GAN=0.5*(F1(I,J+1,K)+F1(I,J,K))
      GAS=0.5*(F1(I,J-1,K)+F1(I,J,K))
      GAT=0.5*(F1(I,J,K+1)+F1(I,J,K))
      GAB=0.5*(F1(I,J,K-1)+F1(I,J,K))
      DENE=0.5*(DEN(I+1,J,K)+DEN(I,J,K))
      DENW=0.5*(DEN(I-1,J,K)+DEN(I,J,K))
      DENN=0.5*(DEN(I,J+1,K)+DEN(I,J,K))
      DENS=0.5*(DEN(I,J-1,K)+DEN(I,J,K))
      DENT=0.5*(DEN(I,J,K+1)+DEN(I,J,K))
      DENB=0.5*(DEN(I,J,K-1)+DEN(I,J,K))
      UE=0.5*(U(I+1,J,K)+U(I,J,K))
      UW=0.5*(U(I-1,J,K)+U(I,J,K))
      UN=0.5*(U(I,J+1,K)+U(I,J,K))
      US=0.5*(U(I,J-1,K)+U(I,J,K))
      UT=0.5*(U(I,J,K+1)+U(I,J,K))
      UB=0.5*(U(I,J,K-1)+U(I,J,K))
      VE=0.5*(V(I+1,J,K)+V(I,J,K))
      VW=0.5*(V(I-1,J,K)+V(I,J,K))
      VN=0.5*(V(I,J+1,K)+V(I,J,K))
      VS=0.5*(V(I,J-1,K)+V(I,J,K))
      VT=0.5*(V(I,J,K+1)+V(I,J,K))
      VB=0.5*(V(I,J,K-1)+V(I,J,K))
      WE=0.5*(W(I+1,J,K)+W(I,J,K))
      WW=0.5*(W(I-1,J,K)+W(I,J,K))
      WN=0.5*(W(I,J+1,K)+W(I,J,K))
      WS=0.5*(W(I,J-1,K)+W(I,J,K))
      WT=0.5*(W(I,J,K+1)+W(I,J,K))
      WB=0.5*(W(I,J,K-1)+W(I,J,K))
      CE=DENE*(UE*CX(I,J,K)+VE*CY(I,J,K)+WE*CZ(I,J,K))
      CW=DENW*(UW*CX(I,J,K)+VW*CY(I,J,K)+WW*CZ(I,J,K))
      CN=DENN*(UN*EX(I,J,K)+VN*EY(I,J,K)+WN*EZ(I,J,K))
      CS=DENS*(US*EX(I,J,K)+VS*EY(I,J,K)+WS*EZ(I,J,K))
      CT=DENT*(UT*SX(I,J,K)+VT*SY(I,J,K)+WT*SZ(I,J,K))
      CB=DENB*(UB*SX(I,J,K)+VB*SY(I,J,K)+WB*SZ(I,J,K))
      DDE=GAE*TXXE(I,J,K)+GAN*TXYN(I,J,K)-GAS*TXYS(I,J,K)+
     1    GAT*TXZT(I,J,K)-GAB*TXZB(I,J,K)
      DDW=GAW*TXXW(I,J,K)+GAS*TXYS(I,J,K)-GAN*TXYN(I,J,K)-
     1    GAT*TXZT(I,J,K)+GAB*TXZB(I,J,K)
      DDN=GAN*TYYN(I,J,K)+GAE*TYXE(I,J,K)-GAW*TYXW(I,J,K)+
     1    GAT*TYZT(I,J,K)-GAB*TYZB(I,J,K)
      DDS=GAS*TYYS(I,J,K)+GAW*TYXW(I,J,K)-GAE*TYXE(I,J,K)-
     1    GAT*TYZT(I,J,K)+GAB*TYZB(I,J,K)
      DDT=GAT*TZZT(I,J,K)+GAE*TZXE(I,J,K)-GAW*TZXW(I,J,K)+
     1    GAN*TZYN(I,J,K)-GAS*TZYS(I,J,K)
      DDB=GAB*TZZB(I,J,K)+GAW*TZXW(I,J,K)-GAE*TZXE(I,J,K)-
     1    GAN*TZYN(I,J,K)+GAS*TZYS(I,J,K)
      CPO=ABS(CE-CW+CN-CS+CT-CB)
      AE(I,J,K)=(AMAX1(ABS(0.5*CE),DDE)-0.5*CE)*TJO(I,J,K)
      AW(I,J,K)=(AMAX1(ABS(0.5*CW),DDW)+0.5*CW)*TJO(I,J,K)
      AN(I,J,K)=(AMAX1(ABS(0.5*CN),DDN)-0.5*CN)*TJO(I,J,K)
      AS(I,J,K)=(AMAX1(ABS(0.5*CS),DDS)+0.5*CS)*TJO(I,J,K)
      AT(I,J,K)=(AMAX1(ABS(0.5*CT),DDT)-0.5*CT)*TJO(I,J,K)
      AB(I,J,K)=(AMAX1(ABS(0.5*CB),DDB)+0.5*CB)*TJO(I,J,K)
      DDNE=GAE*TYXE(I,J,K)+GAN*TXYN(I,J,K)
      DDSE=-GAE*TYXE(I,J,K)-GAS*TXYS(I,J,K)
      DDNW=-GAW*TYXW(I,J,K)-GAN*TXYN(I,J,K)
      DDSW=GAW*TYXW(I,J,K)+GAS*TXYS(I,J,K)
      DDET=GAE*TZXE(I,J,K)+GAT*TXZT(I,J,K)
      DDEB=-GAE*TZXE(I,J,K)-GAB*TXZB(I,J,K)
      DDWT=-GAW*TZXW(I,J,K)-GAT*TXZT(I,J,K)
      DDWB=GAW*TZXW(I,J,K)+GAB*TXZB(I,J,K)
      DDNT=GAN*TZYN(I,J,K)+GAT*TYZT(I,J,K)
      DDNB=-GAN*TZYN(I,J,K)-GAB*TYZB(I,J,K)
      DDST=-GAS*TZYS(I,J,K)-GAT*TYZT(I,J,K)
      DDSB=GAS*TZYS(I,J,K)+GAB*TYZB(I,J,K)
      SU(I,J,K)=CPO*F(I,J,K)+DDNE*F(I+1,J+1,K)+DDSE*F(I+1,J-1,K)+
     1         DDNW*F(I-1,J+1,K)+DDSW*F(I-1,J-1,K)+
     2         DDET*F(I+1,J,K+1)+DDEB*F(I+1,J,K-1)+
     3         DDWT*F(I-1,J,K+1)+DDWB*F(I-1,J,K-1)+
     4         DDNT*F(I,J+1,K+1)+DDNB*F(I,J+1,K-1)+
     5         DDST*F(I,J-1,K+1)+DDSB*F(I,J-1,K-1)
      DV(I,J,K)=CPO*V(I,J,K)+DDNE*V(I+1,J+1,K)+DDSE*V(I+1,J-1,K)+
     1         DDNW*V(I-1,J+1,K)+DDSW*V(I-1,J-1,K)+
     2         DDET*V(I+1,J,K+1)+DDEB*V(I+1,J,K-1)+
     3         DDWT*V(I-1,J,K+1)+DDWB*V(I-1,J,K-1)+
     4         DDNT*V(I,J+1,K+1)+DDNB*V(I,J+1,K-1)+
     5         DDST*V(I,J-1,K+1)+DDSB*V(I,J-1,K-1)
      DW(I,J,K)=CPO*W(I,J,K)+DDNE*W(I+1,J+1,K)+DDSE*W(I+1,J-1,K)+
     1         DDNW*W(I-1,J+1,K)+DDSW*W(I-1,J-1,K)+
     2         DDET*W(I+1,J,K+1)+DDEB*W(I+1,J,K-1)+
     3         DDWT*W(I-1,J,K+1)+DDWB*W(I-1,J,K-1)+
     4         DDNT*W(I,J+1,K+1)+DDNB*W(I,J+1,K-1)+
     5         DDST*W(I,J-1,K+1)+DDSB*W(I,J-1,K-1)
      PP(I,J,K)=CPO*DK(I,J,K)+DDNE*DK(I+1,J+1,K)+DDSE*DK(I+1,J-1,K)+
     1         DDNW*DK(I-1,J+1,K)+DDSW*DK(I-1,J-1,K)+
     2         DDET*DK(I+1,J,K+1)+DDEB*DK(I+1,J,K-1)+
     3         DDWT*DK(I-1,J,K+1)+DDWB*DK(I-1,J,K-1)+
     4         DDNT*DK(I,J+1,K+1)+DDNB*DK(I,J+1,K-1)+
     5         DDST*DK(I,J-1,K+1)+DDSB*DK(I,J-1,K-1)
      SUK(I,J,K)=CPO
      SPK(I,J,K)=-CPO-APO(I,J,K)
      UCXI=UE-UW
      UEDA=UN-US
      USCI=UT-UB
      VCXI=VE-VW
      VEDA=VN-VS
      VSCI=VT-VB
      WCXI=WE-WW
      WEDA=WN-WS
      WSCI=WT-WB
      UX=UCXI*CX(I,J,K)+UEDA*EX(I,J,K)+USCI*SX(I,J,K)
      UY=UCXI*CY(I,J,K)+UEDA*EY(I,J,K)+USCI*SY(I,J,K)
      UZ=UCXI*CZ(I,J,K)+UEDA*EZ(I,J,K)+USCI*SZ(I,J,K)
      VX=VCXI*CX(I,J,K)+VEDA*EX(I,J,K)+VSCI*SX(I,J,K)
      VY=VCXI*CY(I,J,K)+VEDA*EY(I,J,K)+VSCI*SY(I,J,K)
      VZ=VCXI*CZ(I,J,K)+VEDA*EZ(I,J,K)+VSCI*SZ(I,J,K)
      WX=WCXI*CX(I,J,K)+WEDA*EX(I,J,K)+WSCI*SX(I,J,K)
      WY=WCXI*CY(I,J,K)+WEDA*EY(I,J,K)+WSCI*SY(I,J,K)
      WZ=WCXI*CZ(I,J,K)+WEDA*EZ(I,J,K)+WSCI*SZ(I,J,K)
      GEN(I,J,K)=VISE(I,J,K)*((UY+VX)**2+(VZ+WY)**2+(WX+UZ)**2+
     1         2*(UX*UX+VY*VY+WZ*WZ))
  20  CONTINUE
C---------------------------------------------------------------------
  29  CONTINUE
C-----CALCULATE SOURCE TERMS
      GO TO (1,2,3,4,5,6,7), IE
C-----U-, V-, W-SOURCES
  1   CONTINUE
      DO 15 I=IS,IT
      DO 15 J=JS,JT
      DO 15 K=KS,KT
      GAE=0.5*(F1(I+1,J,K)+F1(I,J,K))
      GAW=0.5*(F1(I-1,J,K)+F1(I,J,K))
      GAN=0.5*(F1(I,J+1,K)+F1(I,J,K))
      GAS=0.5*(F1(I,J-1,K)+F1(I,J,K))
      GAT=0.5*(F1(I,J,K+1)+F1(I,J,K))
      GAB=0.5*(F1(I,J,K-1)+F1(I,J,K))
      PCXI=(P(I+1,J,K)+P(I+1,J+1,K)+P(I+1,J,K+1)+P(I+1,J+1,K+1)-
     1      P(I,J,K)-P(I,J+1,K)-P(I,J,K+1)-P(I,J+1,K+1))*0.25
      PEDA=(P(I,J+1,K)+P(I+1,J+1,K)+P(I,J+1,K+1)+P(I+1,J+1,K+1)-
     1      P(I,J,K)-P(I+1,J,K)-P(I,J,K+1)-P(I+1,J,K+1))*0.25
      PSCI=(P(I,J,K+1)+P(I+1,J,K+1)+P(I,J+1,K+1)+P(I+1,J+1,K+1)-
     1      P(I,J,K)-P(I+1,J,K)-P(I,J+1,K)-P(I+1,J+1,K))*0.25
      UCE=U(I+1,J,K)-U(I,J,K)
      UCW=U(I,J,K)-U(I-1,J,K)
      UCN=(U(I+1,J+1,K)+U(I+1,J,K)-U(I-1,J+1,K)-U(I-1,J,K))*0.25
      UCS=(U(I+1,J,K)+U(I+1,J-1,K)-U(I-1,J,K)-U(I-1,J-1,K))*0.25
      UCT=(U(I+1,J,K+1)+U(I+1,J,K)-U(I-1,J,K+1)-U(I-1,J,K))*0.25
      UCB=(U(I+1,J,K)+U(I+1,J,K-1)-U(I-1,J,K)-U(I-1,J,K-1))*0.25
      UEE=(U(I+1,J+1,K)+U(I,J+1,K)-U(I+1,J-1,K)-U(I,J-1,K))*0.25
      UEW=(U(I,J+1,K)+U(I-1,J+1,K)-U(I,J-1,K)-U(I-1,J-1,K))*0.25
      UEN=U(I,J+1,K)-U(I,J,K)
      UES=U(I,J,K)-U(I,J-1,K)
      UET=(U(I,J+1,K+1)+U(I,J+1,K)-U(I,J-1,K+1)-U(I,J-1,K))*0.25
      UEB=(U(I,J+1,K)+U(I,J+1,K-1)-U(I,J-1,K)-U(I,J-1,K-1))*0.25
      USE=(U(I+1,J,K+1)+U(I,J,K+1)-U(I+1,J,K-1)-U(I,J,K-1))*0.25
      USW=(U(I,J,K+1)+U(I-1,J,K+1)-U(I,J,K-1)-U(I-1,J,K-1))*0.25
      USN=(U(I,J+1,K+1)+U(I,J,K+1)-U(I,J+1,K-1)-U(I,J,K-1))*0.25
      USS=(U(I,J,K+1)+U(I,J-1,K+1)-U(I,J,K-1)-U(I,J-1,K-1))*0.25
      UST=U(I,J,K+1)-U(I,J,K)
      USB=U(I,J,K)-U(I,J,K-1)
      VCE=V(I+1,J,K)-V(I,J,K)
      VCW=V(I,J,K)-V(I-1,J,K)
      VCN=(V(I+1,J+1,K)+V(I+1,J,K)-V(I-1,J+1,K)-V(I-1,J,K))*0.25
      VCS=(V(I+1,J,K)+V(I+1,J-1,K)-V(I-1,J,K)-V(I-1,J-1,K))*0.25
      VCT=(V(I+1,J,K+1)+V(I+1,J,K)-V(I-1,J,K+1)-V(I-1,J,K))*0.25
      VCB=(V(I+1,J,K)+V(I+1,J,K-1)-V(I-1,J,K)-V(I-1,J,K-1))*0.25
      VEE=(V(I+1,J+1,K)+V(I,J+1,K)-V(I+1,J-1,K)-V(I,J-1,K))*0.25
      VEW=(V(I,J+1,K)+V(I-1,J+1,K)-V(I,J-1,K)-V(I-1,J-1,K))*0.25
      VEN=V(I,J+1,K)-V(I,J,K)
      VES=V(I,J,K)-V(I,J-1,K)
      VET=(V(I,J+1,K+1)+V(I,J+1,K)-V(I,J-1,K+1)-V(I,J-1,K))*0.25
      VEB=(V(I,J+1,K)+V(I,J+1,K-1)-V(I,J-1,K)-V(I,J-1,K-1))*0.25
      VSE=(V(I+1,J,K+1)+V(I,J,K+1)-V(I+1,J,K-1)-V(I,J,K-1))*0.25
      VSW=(V(I,J,K+1)+V(I-1,J,K+1)-V(I,J,K-1)-V(I-1,J,K-1))*0.25
      VSN=(V(I,J+1,K+1)+V(I,J,K+1)-V(I,J+1,K-1)-V(I,J,K-1))*0.25
      VSS=(V(I,J,K+1)+V(I,J-1,K+1)-V(I,J,K-1)-V(I,J-1,K-1))*0.25
      VST=V(I,J,K+1)-V(I,J,K)
      VSB=V(I,J,K)-V(I,J,K-1)
      WCE=W(I+1,J,K)-W(I,J,K)
      WCW=W(I,J,K)-W(I-1,J,K)
      WCN=(W(I+1,J+1,K)+W(I+1,J,K)-W(I-1,J+1,K)-W(I-1,J,K))*0.25
      WCS=(W(I+1,J,K)+W(I+1,J-1,K)-W(I-1,J,K)-W(I-1,J-1,K))*0.25
      WCT=(W(I+1,J,K+1)+W(I+1,J,K)-W(I-1,J,K+1)-W(I-1,J,K))*0.25
      WCB=(W(I+1,J,K)+W(I+1,J,K-1)-W(I-1,J,K)-W(I-1,J,K-1))*0.25
      WEE=(W(I+1,J+1,K)+W(I,J+1,K)-W(I+1,J-1,K)-W(I,J-1,K))*0.25
      WEW=(W(I,J+1,K)+W(I-1,J+1,K)-W(I,J-1,K)-W(I-1,J-1,K))*0.25
      WEN=W(I,J+1,K)-W(I,J,K)
      WES=W(I,J,K)-W(I,J-1,K)
      WET=(W(I,J+1,K+1)+W(I,J+1,K)-W(I,J-1,K+1)-W(I,J-1,K))*0.25
      WEB=(W(I,J+1,K)+W(I,J+1,K-1)-W(I,J-1,K)-W(I,J-1,K-1))*0.25
      WSE=(W(I+1,J,K+1)+W(I,J,K+1)-W(I+1,J,K-1)-W(I,J,K-1))*0.25
      WSW=(W(I,J,K+1)+W(I-1,J,K+1)-W(I,J,K-1)-W(I-1,J,K-1))*0.25
      WSN=(W(I,J+1,K+1)+W(I,J,K+1)-W(I,J+1,K-1)-W(I,J,K-1))*0.25
      WSS=(W(I,J,K+1)+W(I,J-1,K+1)-W(I,J,K-1)-W(I,J-1,K-1))*0.25
      WST=W(I,J,K+1)-W(I,J,K)
      WSB=W(I,J,K)-W(I,J,K-1)
      SU(I,J,K)=SU(I,J,K)-PCXI*CX(I,J,K)-PEDA*EX(I,J,K)-PSCI*SX(I,J,K)
      DV(I,J,K)=DV(I,J,K)-PCXI*CY(I,J,K)-PEDA*EY(I,J,K)-PSCI*SY(I,J,K)
      DW(I,J,K)=DW(I,J,K)-PCXI*CZ(I,J,K)-PEDA*EZ(I,J,K)-PSCI*SZ(I,J,K)
      CXE=(CX(I+1,J,K)+CX(I,J,K))*0.5
      CXW=(CX(I-1,J,K)+CX(I,J,K))*0.5
      CXN=(CX(I,J+1,K)+CX(I,J,K))*0.5
      CXS=(CX(I,J,K)+CX(I,J-1,K))*0.5
      CXT=(CX(I,J,K+1)+CX(I,J,K))*0.5
      CXB=(CX(I,J,K)+CX(I,J,K-1))*0.5
      EXE=(EX(I+1,J,K)+EX(I,J,K))*0.5
      EXW=(EX(I,J,K)+EX(I-1,J,K))*0.5
      EXN=(EX(I,J+1,K)+EX(I,J,K))*0.5
      EXS=(EX(I,J,K)+EX(I,J-1,K))*0.5
      EXT=(EX(I,J,K+1)+EX(I,J,K))*0.5
      EXB=(EX(I,J,K)+EX(I,J,K-1))*0.5
      SXE=(SX(I+1,J,K)+SX(I,J,K))*0.5
      SXW=(SX(I,J,K)+SX(I-1,J,K))*0.5
      SXN=(SX(I,J+1,K)+SX(I,J,K))*0.5
      SXS=(SX(I,J,K)+SX(I,J-1,K))*0.5
      SXT=(SX(I,J,K+1)+SX(I,J,K))*0.5
      SXB=(SX(I,J,K)+SX(I,J,K-1))*0.5
      QE=GAE*(UCE*CXE+UEE*EXE+USE*SXE)
      QW=GAW*(UCW*CXW+UEW*EXW+USW*SXW)
      QN=GAN*(UCN*CXN+UEN*EXN+USN*SXN)
      QS=GAS*(UCS*CXS+UES*EXS+USS*SXS)
      QT=GAT*(UCT*CXT+UET*EXT+UST*SXT)
      QB=GAB*(UCB*CXB+UEB*EXB+USB*SXB)
      SOC1=CX(I,J,K)*(QE-QW)+EX(I,J,K)*(QN-QS)+SX(I,J,K)*(QT-QB)
      QE=GAE*(VCE*CXE+VEE*EXE+VSE*SXE)
      QW=GAW*(VCW*CXW+VEW*EXW+VSW*SXW)
      QN=GAN*(VCN*CXN+VEN*EXN+VSN*SXN)
      QS=GAS*(VCS*CXS+VES*EXS+VSS*SXS)
      QT=GAT*(VCT*CXT+VET*EXT+VST*SXT)
      QB=GAB*(VCB*CXB+VEB*EXB+VSB*SXB)
      SOC2=CY(I,J,K)*(QE-QW)+EY(I,J,K)*(QN-QS)+SY(I,J,K)*(QT-QB)
      QE=GAE*(WCE*CXE+WEE*EXE+WSE*SXE)
      QW=GAW*(WCW*CXW+WEW*EXW+WSW*SXW)
      QN=GAN*(WCN*CXN+WEN*EXN+WSN*SXN)
      QS=GAS*(WCS*CXS+WES*EXS+WSS*SXS)
      QT=GAT*(WCT*CXT+WET*EXT+WST*SXT)
      QB=GAB*(WCB*CXB+WEB*EXB+WSB*SXB)
      SOC3=CZ(I,J,K)*(QE-QW)+EZ(I,J,K)*(QN-QS)+SZ(I,J,K)*(QT-QB)
      SU(I,J,K)=SU(I,J,K)+SOC1+SOC2+SOC3+APO(I,J,K)*FO(I,J,K)
      CYE=(CY(I+1,J,K)+CY(I,J,K))*0.5
      CYW=(CY(I,J,K)+CY(I-1,J,K))*0.5
      CYN=(CY(I,J+1,K)+CY(I,J,K))*0.5
      CYS=(CY(I,J,K)+CY(I,J-1,K))*0.5
      CYT=(CY(I,J,K+1)+CY(I,J,K))*0.5
      CYB=(CY(I,J,K)+CY(I,J,K-1))*0.5
      EYE=(EY(I+1,J,K)+EY(I,J,K))*0.5
      EYW=(EY(I,J,K)+EY(I-1,J,K))*0.5
      EYN=(EY(I,J+1,K)+EY(I,J,K))*0.5
      EYS=(EY(I,J,K)+EY(I,J-1,K))*0.5
      EYT=(EY(I,J,K+1)+EY(I,J,K))*0.5
      EYB=(EY(I,J,K)+EY(I,J,K-1))*0.5
      SYE=(SY(I+1,J,K)+SY(I,J,K))*0.5
      SYW=(SY(I,J,K)+SY(I-1,J,K))*0.5
      SYN=(SY(I,J+1,K)+SY(I,J,K))*0.5
      SYS=(SY(I,J,K)+SY(I,J-1,K))*0.5
      SYT=(SY(I,J,K+1)+SY(I,J,K))*0.5
      SYB=(SY(I,J,K)+SY(I,J,K-1))*0.5
      QE=GAE*(UCE*CYE+UEE*EYE+USE*SYE)
      QW=GAW*(UCW*CYW+UEW*EYW+USW*SYW)
      QN=GAN*(UCN*CYN+UEN*EYN+USN*SYN)
      QS=GAS*(UCS*CYS+UES*EYS+USS*SYS)
      QT=GAT*(UCT*CYT+UET*EYT+UST*SYT)
      QB=GAB*(UCB*CYB+UEB*EYB+USB*SYB)
      SOC1=CX(I,J,K)*(QE-QW)+EX(I,J,K)*(QN-QS)+SX(I,J,K)*(QT-QB)
      QE=GAE*(VCE*CYE+VEE*EYE+VSE*SYE)
      QW=GAW*(VCW*CYW+VEW*EYW+VSW*SYW)
      QN=GAN*(VCN*CYN+VEN*EYN+VSN*SYN)
      QS=GAS*(VCS*CYS+VES*EYS+VSS*SYS)
      QT=GAT*(VCT*CYT+VET*EYT+VST*SYT)
      QB=GAB*(VCB*CYB+VEB*EYB+VSB*SYB)
      SOC2=CY(I,J,K)*(QE-QW)+EY(I,J,K)*(QN-QS)+SY(I,J,K)*(QT-QB)
      QE=GAE*(WCE*CYE+WEE*EYE+WSE*SYE)
      QW=GAW*(WCW*CYW+WEW*EYW+WSW*SYW)
      QN=GAN*(WCN*CYN+WEN*EYN+WSN*SYN)
      QS=GAS*(WCS*CYS+WES*EYS+WSS*SYS)
      QT=GAT*(WCT*CYT+WET*EYT+WST*SYT)
      QB=GAB*(WCB*CYB+WEB*EYB+WSB*SYB)
      SOC3=CZ(I,J,K)*(QE-QW)+EZ(I,J,K)*(QN-QS)+SZ(I,J,K)*(QT-QB)
      DV(I,J,K)=DV(I,J,K)+SOC1+SOC2+SOC3
      CZE=(CZ(I+1,J,K)+CZ(I,J,K))*0.5
      CZW=(CZ(I,J,K)+CZ(I-1,J,K))*0.5
      CZN=(CZ(I,J+1,K)+CZ(I,J,K))*0.5
      CZS=(CZ(I,J,K)+CZ(I,J-1,K))*0.5
      CZT=(CZ(I,J,K+1)+CZ(I,J,K))*0.5
      CZB=(CZ(I,J,K)+CZ(I,J,K-1))*0.5
      EZE=(EZ(I+1,J,K)+EZ(I,J,K))*0.5
      EZW=(EZ(I,J,K)+EZ(I-1,J,K))*0.5
      EZN=(EZ(I,J+1,K)+EZ(I,J,K))*0.5
      EZS=(EZ(I,J,K)+EZ(I,J-1,K))*0.5
      EZT=(EZ(I,J,K+1)+EZ(I,J,K))*0.5
      EZB=(EZ(I,J,K)+EZ(I,J,K-1))*0.5
      SZE=(SZ(I+1,J,K)+SZ(I,J,K))*0.5
      SZW=(SZ(I,J,K)+SZ(I-1,J,K))*0.5
      SZN=(SZ(I,J+1,K)+SZ(I,J,K))*0.5
      SZS=(SZ(I,J,K)+SZ(I,J-1,K))*0.5
      SZT=(SZ(I,J,K+1)+SZ(I,J,K))*0.5
      SZB=(SZ(I,J,K)+SZ(I,J,K-1))*0.5
      QE=GAE*(UCE*CZE+UEE*EZE+USE*SZE)
      QW=GAW*(UCW*CZW+UEW*EZW+USW*SZW)
      QN=GAN*(UCN*CZN+UEN*EZN+USN*SZN)
      QS=GAS*(UCS*CZS+UES*EZS+USS*SZS)
      QT=GAT*(UCT*CZT+UET*EZT+UST*SZT)
      QB=GAB*(UCB*CZB+UEB*EZB+USB*SZB)
      SOC1=CX(I,J,K)*(QE-QW)+EX(I,J,K)*(QN-QS)+SX(I,J,K)*(QT-QB)
      QE=GAE*(VCE*CZE+VEE*EZE+VSE*SZE)
      QW=GAW*(VCW*CZW+VEW*EZW+VSW*SZW)
      QN=GAN*(VCN*CZN+VEN*EZN+VSN*SZN)
      QS=GAS*(VCS*CZS+VES*EZS+VSS*SZS)
      QT=GAT*(VCT*CZT+VET*EZT+VST*SZT)
      QB=GAB*(VCB*CZB+VEB*EZB+VSB*SZB)
      SOC2=CY(I,J,K)*(QE-QW)+EY(I,J,K)*(QN-QS)+SY(I,J,K)*(QT-QB)
      QE=GAE*(WCE*CZE+WEE*EZE+WSE*SZE)
      QW=GAW*(WCW*CZW+WEW*EZW+WSW*SZW)
      QN=GAN*(WCN*CZN+WEN*EZN+WSN*SZN)
      QS=GAS*(WCS*CZS+WES*EZS+WSS*SZS)
      QT=GAT*(WCT*CZT+WET*EZT+WST*SZT)
      QB=GAB*(WCB*CZB+WEB*EZB+WSB*SZB)
      SOC3=CZ(I,J,K)*(QE-QW)+EZ(I,J,K)*(QN-QS)+SZ(I,J,K)*(QT-QB)
      DW(I,J,K)=DW(I,J,K)+SOC1+SOC2+SOC3
      SU(I,J,K)=SU(I,J,K)*TJO(I,J,K)
      SP(I,J,K)=SPK(I,J,K)*TJO(I,J,K)
  15  CONTINUE
      GO TO 60
C-----V-SOURCE
  2   CONTINUE
      DO 25 I=IS,IT
      DO 25 J=JS,JT
      DO 25 K=KS,KT
      SU(I,J,K)=(DV(I,J,K)+APO(I,J,K)*FO(I,J,K))*TJO(I,J,K)
  25  CONTINUE
      GO TO 60
C-----W-SOURCE
  3   CONTINUE
      DO 35 I=IS,IT
      DO 35 J=JS,JT
      DO 35 K=KS,KT
      SU(I,J,K)=(DW(I,J,K)+APO(I,J,K)*FO(I,J,K))*TJO(I,J,K)
  35  CONTINUE
      GO TO 60
C-----TM-SOURCE
  4   CONTINUE
      DO 45 I=IS,IT
      DO 45 J=JS,JT
      DO 45 K=KS,KT
      GAE=0.5*(VISE(I+1,J,K)+VISE(I,J,K))/SIGF
      GAW=0.5*(VISE(I-1,J,K)+VISE(I,J,K))/SIGF
      GAN=0.5*(VISE(I,J+1,K)+VISE(I,J,K))/SIGF
      GAS=0.5*(VISE(I,J-1,K)+VISE(I,J,K))/SIGF
      GAT=0.5*(VISE(I,J,K+1)+VISE(I,J,K))/SIGF
      GAB=0.5*(VISE(I,J,K-1)+VISE(I,J,K))/SIGF
      DDNE=GAE*TYXE(I,J,K)+GAN*TXYN(I,J,K)
      DDSE=-GAE*TYXE(I,J,K)-GAS*TXYS(I,J,K)
      DDNW=-GAW*TYXW(I,J,K)-GAN*TXYN(I,J,K)
      DDSW=GAW*TYXW(I,J,K)+GAS*TXYS(I,J,K)
      DDET=GAE*TZXE(I,J,K)+GAT*TXZT(I,J,K)
      DDEB=-GAE*TZXE(I,J,K)-GAB*TXZB(I,J,K)
      DDWT=-GAW*TZXW(I,J,K)-GAT*TXZT(I,J,K)
      DDWB=GAW*TZXW(I,J,K)+GAB*TXZB(I,J,K)
      DDNT=GAN*TZYN(I,J,K)+GAT*TYZT(I,J,K)
      DDNB=-GAN*TZYN(I,J,K)-GAB*TYZB(I,J,K)
      DDST=-GAS*TZYS(I,J,K)-GAT*TYZT(I,J,K)
      DDSB=GAS*TZYS(I,J,K)+GAB*TYZB(I,J,K)
      SU(I,J,K)=DDNE*F(I+1,J+1,K)+DDSE*F(I+1,J-1,K)+
     1         DDNW*F(I-1,J+1,K)+DDSW*F(I-1,J-1,K)+
     2         DDET*F(I+1,J,K+1)+DDEB*F(I+1,J,K-1)+
     3         DDWT*F(I-1,J,K+1)+DDWB*F(I-1,J,K-1)+
     4         DDNT*F(I,J+1,K+1)+DDNB*F(I,J+1,K-1)+
     5         DDST*F(I,J-1,K+1)+DDSB*F(I,J-1,K-1)
      SU(I,J,K)=(SU(I,J,K)+APO(I,J,K)*FO(I,J,K))*TJO(I,J,K)
      SP(I,J,K)=-APO(I,J,K)*TJO(I,J,K)
  45  CONTINUE
      GO TO 60
C-----K-SOURCE
  5   CONTINUE
      DO 55 I=IS,IT
      DO 55 J=JS,JT
      DO 55 K=KS,KT
      SU(I,J,K)=GEN(I,J,K)+PP(I,J,K)+APO(I,J,K)*FO(I,J,K)
      P1=DEN(I,J,K)**2
      SP(I,J,K)=SPK(I,J,K)-CMU*P1*F(I,J,K)/VISE(I,J,K)
      SU(I,J,K)=SU(I,J,K)*TJO(I,J,K)
      SP(I,J,K)=SP(I,J,K)*TJO(I,J,K)
  55  CONTINUE
      GO TO 60
C-----E-SOURCE
  6   CONTINUE
      DO 65 I=IS,IT
      DO 65 J=JS,JT
      DO 65 K=KS,KT
      P1=DEN(I,J,K)**2
      SU(I,J,K)=SU(I,J,K)+C1*CMU*GEN(I,J,K)*P1*DK(I,J,K)/
     1          VISE(I,J,K)+APO(I,J,K)*FO(I,J,K)
      TMDK=DK(I,J,K)+SMNUM
      SP(I,J,K)=SPK(I,J,K)-C2*DEN(I,J,K)*F(I,J,K)/TMDK
      SU(I,J,K)=SU(I,J,K)*TJO(I,J,K)
      SP(I,J,K)=SP(I,J,K)*TJO(I,J,K)
  65  CONTINUE
      GO TO 60
C-----PF-SOURCE
  7   CONTINUE
      DO 75 I=IS,IT
      DO 75 J=JS,JT
      DO 75 K=KS,KT
      SU(I,J,K)=SU(I,J,K)+0.85*GEN(I,J,K)+APO(I,J,K)*FO(I,J,K)
      TMDK=DK(I,J,K)+SMNUM
      SP(I,J,K)=SPK(I,J,K)-DE(I,J,K)/TMDK
      SU(I,J,K)=SU(I,J,K)*TJO(I,J,K)
      SP(I,J,K)=SP(I,J,K)*TJO(I,J,K)
  75  CONTINUE
  60  CONTINUE
C-----MODIFY WALL BOUNDARY CONDITIONS THRU WALL FUNCTIONS
      IF(IG .NE. 2) GO TO 410
      CALL BOUNC(IE,F)
 410  CONTINUE
C-----SET SYMMETRIC, CYCLIC AND EXIT LINK COEFF.
      CALL SYMOUT(2,IE,IS,IT,JS,JT,KS,KT)
C-----LINK COEFF. ASSEMBLY AND BLOCKAGES
      DO 500 I=IS,IT
      DO 500 J=JS,JT
      DO 500 K=KS,KT
      F1(I,J,K)=F(I,J,K)
      ANAB=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)+
     1     AB(I,J,K)
      AP(I,J,K)=ANAB-SP(I,J,K)
      PDUV=1.0
      IF(MC(I,J,K) .LT. 1) GO TO 530
      AP(I,J,K)=ALF
      AN(I,J,K)=0.0
      AS(I,J,K)=0.0
      AE(I,J,K)=0.0
      AW(I,J,K)=0.0
      AT(I,J,K)=0.0
      AB(I,J,K)=0.0
      SU(I,J,K)=F(I,J,K)
      PDUV=0.0
 530  CONTINUE
C-----UNDER-RELAXATION
      P1=1.2*AP(I,J,K)
      AP(I,J,K)=AP(I,J,K)/ALF
      SU(I,J,K)=SU(I,J,K)+PDUV*(1.0-ALF)*AP(I,J,K)*F(I,J,K)
      IF(IE .EQ. 1) DU(I,J,K)=TJO(I,J,K)*PDUV/(P1-ANAB)
 500  CONTINUE
C-----LINEAR EQUATIONS SLOVER
      DO 550 I=1,ISWF
 550  CALL LINERX(1,IS,JS,KS,IT,JT,KT,F)
C-----CALCULATE MAXIMUM CORRECTION OF CURRENT ITERATION
      DO 555 I=IS,IT
      DO 555 J=JS,JT
      DO 555 K=KS,KT
      P1=ABS(F(I,J,K)-F1(I,J,K))
      ERRF=AMAX1(ERRF,P1)
 555  CONTINUE
      RETURN
C---------------------------------------------------------------------
C-----P EQUATION
  10  CONTINUE
      IS=2
      IT=L
      JS=2
      JT=M
      KS=2
      KT=N
C-----SOURCE AT PX LOCATIONS
      DO 310 I=IS,IT
      DO 310 J=JS,JT
      DO 310 K=KS,KT
      F(I,J,K)=0.0
      SUK(I,J,K)=0.0
      SPK(I,J,K)=0.0
      SP(I,J,K)=0.0
      DENE=DEN(I,J,K)
      DENW=DEN(I-1,J,K)
      DENN=0.25*(DEN(I,J,K)+DEN(I-1,J,K)+DEN(I,J+1,K)+DEN(I-1,J+1,K))
      DENS=0.25*(DEN(I,J,K)+DEN(I-1,J,K)+DEN(I,J-1,K)+DEN(I-1,J-1,K))
      DENT=0.25*(DEN(I,J,K)+DEN(I-1,J,K)+DEN(I,J,K+1)+DEN(I-1,J,K+1))
      DENB=0.25*(DEN(I,J,K)+DEN(I-1,J,K)+DEN(I,J,K-1)+DEN(I-1,J,K-1))
      UE=U(I,J,K)
      UW=U(I-1,J,K)
      UN=0.25*(U(I,J,K)+U(I-1,J,K)+U(I,J+1,K)+U(I-1,J+1,K))
      US=0.25*(U(I,J,K)+U(I-1,J,K)+U(I,J-1,K)+U(I-1,J-1,K))
      UT=0.25*(U(I,J,K)+U(I-1,J,K)+U(I,J,K+1)+U(I-1,J,K+1))
      UB=0.25*(U(I,J,K)+U(I-1,J,K)+U(I,J,K-1)+U(I-1,J,K-1))
      VE=V(I,J,K)
      VW=V(I-1,J,K)
      VN=0.25*(V(I,J,K)+V(I-1,J,K)+V(I,J+1,K)+V(I-1,J+1,K))
      VS=0.25*(V(I,J,K)+V(I-1,J,K)+V(I,J-1,K)+V(I-1,J-1,K))
      VT=0.25*(V(I,J,K)+V(I-1,J,K)+V(I,J,K+1)+V(I-1,J,K+1))
      VB=0.25*(V(I,J,K)+V(I-1,J,K)+V(I,J,K-1)+V(I-1,J,K-1))
      WE=W(I,J,K)
      WW=W(I-1,J,K)
      WN=0.25*(W(I,J,K)+W(I-1,J,K)+W(I,J+1,K)+W(I-1,J+1,K))
      WS=0.25*(W(I,J,K)+W(I-1,J,K)+W(I,J-1,K)+W(I-1,J-1,K))
      WT=0.25*(W(I,J,K)+W(I-1,J,K)+W(I,J,K+1)+W(I-1,J,K+1))
      WB=0.25*(W(I,J,K)+W(I-1,J,K)+W(I,J,K-1)+W(I-1,J,K-1))
      CXQ=0.5*(CX(I,J,K)+CX(I-1,J,K))
      EXQ=0.5*(EX(I,J,K)+EX(I-1,J,K))
      SXQ=0.5*(SX(I,J,K)+SX(I-1,J,K))
      CYQ=0.5*(CY(I,J,K)+CY(I-1,J,K))
      EYQ=0.5*(EY(I,J,K)+EY(I-1,J,K))
      SYQ=0.5*(SY(I,J,K)+SY(I-1,J,K))
      CZQ=0.5*(CZ(I,J,K)+CZ(I-1,J,K))
      EZQ=0.5*(EZ(I,J,K)+EZ(I-1,J,K))
      SZQ=0.5*(SZ(I,J,K)+SZ(I-1,J,K))
      CE=DENE*(UE*CXQ+VE*CYQ+WE*CZQ)
      CW=DENW*(UW*CXQ+VW*CYQ+WW*CZQ)
      CN=DENN*(UN*EXQ+VN*EYQ+WN*EZQ)
      CS=DENS*(US*EXQ+VS*EYQ+WS*EZQ)
      CT=DENT*(UT*SXQ+VT*SYQ+WT*SZQ)
      CB=DENB*(UB*SXQ+VB*SYQ+WB*SZQ)
      SUK(I,J,K)=-(CE-CW+CN-CS+CT-CB)*(TJO(I-1,J,K)+TJO(I,J,K))*0.5
 310  CONTINUE
C-----SOURCE AT PY LOCATIONS
      DO 311 I=IS,IT
      DO 311 J=JS,JT
      DO 311 K=KS,KT
      DENN=DEN(I,J,K)
      DENS=DEN(I,J-1,K)
      DENE=0.25*(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I+1,J,K)+DEN(I+1,J-1,K))
      DENW=0.25*(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I-1,J,K)+DEN(I-1,J-1,K))
      DENT=0.25*(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I,J,K+1)+DEN(I,J-1,K+1))
      DENB=0.25*(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I,J,K-1)+DEN(I,J-1,K-1))
      UN=U(I,J,K)
      US=U(I,J-1,K)
      UE=0.25*(U(I,J,K)+U(I,J-1,K)+U(I+1,J,K)+U(I+1,J-1,K))
      UW=0.25*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K))
      UT=0.25*(U(I,J,K)+U(I,J-1,K)+U(I,J,K+1)+U(I,J-1,K+1))
      UB=0.25*(U(I,J,K)+U(I,J-1,K)+U(I,J,K-1)+U(I,J-1,K-1))
      VN=V(I,J,K)
      VS=V(I,J-1,K)
      VE=0.25*(V(I,J,K)+V(I,J-1,K)+V(I+1,J,K)+V(I+1,J-1,K))
      VW=0.25*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K))
      VT=0.25*(V(I,J,K)+V(I,J-1,K)+V(I,J,K+1)+V(I,J-1,K+1))
      VB=0.25*(V(I,J,K)+V(I,J-1,K)+V(I,J,K-1)+V(I,J-1,K-1))
      WN=W(I,J,K)
      WS=W(I,J-1,K)
      WE=0.25*(W(I,J,K)+W(I,J-1,K)+W(I+1,J,K)+W(I+1,J-1,K))
      WW=0.25*(W(I,J,K)+W(I,J-1,K)+W(I-1,J,K)+W(I-1,J-1,K))
      WT=0.25*(W(I,J,K)+W(I,J-1,K)+W(I,J,K+1)+W(I,J-1,K+1))
      WB=0.25*(W(I,J,K)+W(I,J-1,K)+W(I,J,K-1)+W(I,J-1,K-1))
      CXQ=0.5*(CX(I,J,K)+CX(I,J-1,K))
      EXQ=0.5*(EX(I,J,K)+EX(I,J-1,K))
      SXQ=0.5*(SX(I,J,K)+SX(I,J-1,K))
      CYQ=0.5*(CY(I,J,K)+CY(I,J-1,K))
      EYQ=0.5*(EY(I,J,K)+EY(I,J-1,K))
      SYQ=0.5*(SY(I,J,K)+SY(I,J-1,K))
      CZQ=0.5*(CZ(I,J,K)+CZ(I,J-1,K))
      EZQ=0.5*(EZ(I,J,K)+EZ(I,J-1,K))
      SZQ=0.5*(SZ(I,J,K)+SZ(I,J-1,K))
      CE=DENE*(UE*CXQ+VE*CYQ+WE*CZQ)
      CW=DENW*(UW*CXQ+VW*CYQ+WW*CZQ)
      CN=DENN*(UN*EXQ+VN*EYQ+WN*EZQ)
      CS=DENS*(US*EXQ+VS*EYQ+WS*EZQ)
      CT=DENT*(UT*SXQ+VT*SYQ+WT*SZQ)
      CB=DENB*(UB*SXQ+VB*SYQ+WB*SZQ)
      SPK(I,J,K)=-(CE-CW+CN-CS+CT-CB)*(TJO(I,J-1,K)+TJO(I,J,K))*0.5
 311  CONTINUE
C-----SOURCE AT PZ LOCATIONS
      DO 312 I=IS,IT
      DO 312 J=JS,JT
      DO 312 K=KS,KT
      DENT=DEN(I,J,K)
      DENB=DEN(I,J,K-1)
      DENN=0.25*(DEN(I,J,K)+DEN(I,J,K-1)+DEN(I,J+1,K)+DEN(I,J+1,K-1))
      DENS=0.25*(DEN(I,J,K)+DEN(I,J,K-1)+DEN(I,J-1,K)+DEN(I,J-1,K-1))
      DENE=0.25*(DEN(I,J,K)+DEN(I,J,K-1)+DEN(I+1,J,K)+DEN(I+1,J,K-1))
      DENW=0.25*(DEN(I,J,K)+DEN(I,J,K-1)+DEN(I-1,J,K)+DEN(I-1,J,K-1))
      UT=U(I,J,K)
      UB=U(I,J,K-1)
      UN=0.25*(U(I,J,K)+U(I,J,K-1)+U(I,J+1,K)+U(I,J+1,K-1))
      US=0.25*(U(I,J,K)+U(I,J,K-1)+U(I,J-1,K)+U(I,J-1,K-1))
      UE=0.25*(U(I,J,K)+U(I,J,K-1)+U(I+1,J,K)+U(I+1,J,K-1))
      UW=0.25*(U(I,J,K)+U(I,J,K-1)+U(I-1,J,K)+U(I-1,J,K-1))
      VT=V(I,J,K)
      VB=V(I,J,K-1)
      VN=0.25*(V(I,J,K)+V(I,J,K-1)+V(I,J+1,K)+V(I,J+1,K-1))
      VS=0.25*(V(I,J,K)+V(I,J,K-1)+V(I,J-1,K)+V(I,J-1,K-1))
      VE=0.25*(V(I,J,K)+V(I,J,K-1)+V(I+1,J,K)+V(I+1,J,K-1))
      VW=0.25*(V(I,J,K)+V(I,J,K-1)+V(I-1,J,K)+V(I-1,J,K-1))
      WT=W(I,J,K)
      WB=W(I,J,K-1)
      WN=0.25*(W(I,J,K)+W(I,J,K-1)+W(I,J+1,K)+W(I,J+1,K-1))
      WS=0.25*(W(I,J,K)+W(I,J,K-1)+W(I,J-1,K)+W(I,J-1,K-1))
      WE=0.25*(W(I,J,K)+W(I,J,K-1)+W(I+1,J,K)+W(I+1,J,K-1))
      WW=0.25*(W(I,J,K)+W(I,J,K-1)+W(I-1,J,K)+W(I-1,J,K-1))
      CXQ=0.5*(CX(I,J,K)+CX(I,J,K-1))
      EXQ=0.5*(EX(I,J,K)+EX(I,J,K-1))
      SXQ=0.5*(SX(I,J,K)+SX(I,J,K-1))
      CYQ=0.5*(CY(I,J,K)+CY(I,J,K-1))
      EYQ=0.5*(EY(I,J,K)+EY(I,J,K-1))
      SYQ=0.5*(SY(I,J,K)+SY(I,J,K-1))
      CZQ=0.5*(CZ(I,J,K)+CZ(I,J,K-1))
      EZQ=0.5*(EZ(I,J,K)+EZ(I,J,K-1))
      SZQ=0.5*(SZ(I,J,K)+SZ(I,J,K-1))
      CE=DENE*(UE*CXQ+VE*CYQ+WE*CZQ)
      CW=DENW*(UW*CXQ+VW*CYQ+WW*CZQ)
      CN=DENN*(UN*EXQ+VN*EYQ+WN*EZQ)
      CS=DENS*(US*EXQ+VS*EYQ+WS*EZQ)
      CT=DENT*(UT*SXQ+VT*SYQ+WT*SZQ)
      CB=DENB*(UB*SXQ+VB*SYQ+WB*SZQ)
      SP(I,J,K)=-(CE-CW+CN-CS+CT-CB)*(TJO(I,J,K-1)+TJO(I,J,K))*0.5
 312  CONTINUE
C-----LINK COEFF. AND SOURCE TERM AT CELL CENTER
      DO 100 I=IS,IT
      DO 100 J=JS,JT
      DO 100 K=KS,KT
      II=I-1
      JJ=J-1
      KK=K-1
      GAE=(DU(I,J,K)+DU(I,JJ,K)+DU(I,J,KK)+DU(I,JJ,KK))*0.25
      GAW=(DU(II,J,K)+DU(II,JJ,K)+DU(II,J,KK)+DU(II,JJ,KK))*0.25
      GAN=(DU(I,J,K)+DU(II,J,K)+DU(I,J,KK)+DU(II,J,KK))*0.25
      GAS=(DU(I,JJ,K)+DU(II,JJ,K)+DU(I,JJ,KK)+DU(II,JJ,KK))*0.25
      GAT=(DU(I,J,K)+DU(II,J,K)+DU(I,JJ,K)+DU(II,JJ,K))*0.25
      GAB=(DU(I,J,KK)+DU(II,J,KK)+DU(I,JJ,KK)+DU(II,JJ,KK))*0.25
      TXXEQ=(TXXE(I,J,K)+TXXE(II,JJ,KK))*0.5
      TXXWQ=(TXXW(I,J,K)+TXXW(II,JJ,KK))*0.5
      TYYNQ=(TYYN(I,J,K)+TYYN(II,JJ,KK))*0.5
      TYYSQ=(TYYS(I,J,K)+TYYS(II,JJ,KK))*0.5
      TZZTQ=(TZZT(I,J,K)+TZZT(II,JJ,KK))*0.5
      TZZBQ=(TZZB(I,J,K)+TZZB(II,JJ,KK))*0.5
      TQ=(TJO(I,J,K)+TJO(II,JJ,KK))*0.5
      AE(I,J,K)=TQ*(GAE*TXXEQ)
      AW(I,J,K)=TQ*(GAW*TXXWQ)
      AN(I,J,K)=TQ*(GAN*TYYNQ)
      AS(I,J,K)=TQ*(GAS*TYYSQ)
      AT(I,J,K)=TQ*(GAT*TZZTQ)
      AB(I,J,K)=TQ*(GAB*TZZBQ)
C-----P' SOURCE
      P1=1.0
      P2=1.0
      P3=1.0
      P4=1.0
      P5=1.0
      P6=1.0
      IF(GAN .EQ. 0.0) P1=0.0
      IF(GAS .EQ. 0.0) P2=0.0
      IF(GAE .EQ. 0.0) P3=0.0
      IF(GAW .EQ. 0.0) P4=0.0
      IF(GAT .EQ. 0.0) P5=0.0
      IF(GAB .EQ. 0.0) P6=0.0
      PQ=0.5*((P1+P2)*(P3+P4+P5+P6)+(P3+P4)*(P1+P2+P5+P6)+
     1        (P5+P6)*(P1+P2+P3+P4))
      IF(PQ .EQ. 0.0) PQ=1.0
      SU(I,J,K)=(P1*P5*SUK(I,J,K)+P5*P2*SUK(I,J-1,K)+
     1           P2*P6*SUK(I,J-1,K-1)+P6*P1*SUK(I,J,K-1)+
     2           P3*P5*SPK(I,J,K)+P5*P4*SPK(I-1,J,K)+
     3           P4*P6*SPK(I-1,J,K-1)+P6*P3*SPK(I,J,K-1)+
     4           P1*P3*SP(I,J,K)+P3*P2*SP(I,J-1,K)+
     5           P2*P4*SP(I-1,J-1,K)+P4*P1*SP(I-1,J,K))/PQ
 100  CONTINUE
C-----LINK COEFF. ASSEMBLY AND BLOCKAGES
      DO 870 I=IS,IT
      DO 870 J=JS,JT
      DO 870 K=KS,KT
      ANAB=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)+AB(I,J,K)
      AP(I,J,K)=ANAB
      IF(MC(I,J,K) .LT. 2) GO TO 870
      AP(I,J,K)=1.0
      AE(I,J,K)=0.0
      AW(I,J,K)=0.0
      AN(I,J,K)=0.0
      AS(I,J,K)=0.0
      AT(I,J,K)=0.0
      AB(I,J,K)=0.0
      SU(I,J,K)=0.0
 870  CONTINUE
C-----LINEAR EQUATIONS SOLVER
      CALL LINERX(2,IS,JS,KS,IT,JT,KT,F)
      DO 850 I=1,ISWF
 850  CALL LINERX(1,IS,JS,KS,IT,JT,KT,F)
C-----PRESSURE AND VELOCITIES CORRECTIONS
      PPREF=F(2,2,2)
      DO 880 I=IS,IT
      DO 880 J=JS,JT
      DO 880 K=KS,KT
      IF(MC(I,J,K) .GE. 2) GO TO 880
      P1=F(I,J,K)-PPREF
      P(I,J,K)=P(I,J,K)+ALP*P1
      ERRF=AMAX1(ERRF,ABS(P1))
 880  CONTINUE
      DO 600 I=2,LT
      DO 600 J=2,MT
      DO 600 K=2,NT
      IF(MC(I,J,K) .GE. 1) GO TO 600
      II=I+1
      JJ=J+1
      KK=K+1
      PE=(F(II,J,K)+F(II,JJ,K)+F(II,J,KK)+F(II,JJ,KK))*0.25
      PW=(F(I,J,K)+F(I,JJ,K)+F(I,J,KK)+F(I,JJ,KK))*0.25
      PN=(F(I,JJ,K)+F(II,JJ,K)+F(I,JJ,KK)+F(II,JJ,KK))*0.25
      P5=(F(I,J,K)+F(II,J,K)+F(I,J,KK)+F(II,J,KK))*0.25
      PT=(F(I,J,KK)+F(II,J,KK)+F(I,JJ,KK)+F(II,JJ,KK))*0.25
      P3=(F(I,J,K)+F(II,J,K)+F(I,JJ,K)+F(II,JJ,K))*0.25
      PCXI=PE-PW
      PEDA=PN-P5
      PSCI=PT-P3
      CXQ=CX(I,J,K)
      EXQ=EX(I,J,K)
      SXQ=SX(I,J,K)
      CYQ=CY(I,J,K)
      EYQ=EY(I,J,K)
      SYQ=SY(I,J,K)
      CZQ=CZ(I,J,K)
      EZQ=EZ(I,J,K)
      SZQ=SZ(I,J,K)
      PXX=PCXI*CXQ+PEDA*EXQ+PSCI*SXQ
      PYY=PCXI*CYQ+PEDA*EYQ+PSCI*SYQ
      PZZ=PCXI*CZQ+PEDA*EZQ+PSCI*SZQ
      U(I,J,K)=U(I,J,K)-DU(I,J,K)*PXX
      V(I,J,K)=V(I,J,K)-DU(I,J,K)*PYY
      W(I,J,K)=W(I,J,K)-DU(I,J,K)*PZZ
 600  CONTINUE
      RETURN
      END
      SUBROUTINE LINERX(ISOL,IS,JS,KS,LT,MT,NT,F)
      DIMENSION A(20),B(20),C(20),D(20),F(21,18,10),BB(21,18,10)
      COMMON
     1/COEF/ AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2   SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3   AS(21,18,10),AT(21,18,10),AB(21,18,10),APO(21,18,10)
C-----LINEAR EQUATIONS SOLVERS
      GO TO (1,2),ISOL
C-----LINE-RELAXATION USING TOMA
  1   CONTINUE
      DO 10 I=IS,LT
      DO 10 J=JS,MT
      DO 10 K=KS,NT
      BB(I,J,K)=AE(I,J,K)*F(I+1,J,K)+AN(I,J,K)*F(I,J+1,K)+
     1          SU(I,J,K)
  10  CONTINUE
      A(KS-1)=0.0
      DO 100 I=IS,LT
      DO 100 J=JS,MT
      C(KS-1)=F(I,J,KS-1)
      DO 101 K=KS,NT
      A(K)=AT(I,J,K)
      B(K)=AB(I,J,K)
      C(K)=BB(I,J,K)+AW(I,J,K)*F(I-1,J,K)+
     1     AS(I,J,K)*F(I,J-1,K)
      D(K)=AP(I,J,K)
      TERM1=D(K)-B(K)*A(K-1)+1.E-30
      TERM=1.0/TERM1
      A(K)=A(K)*TERM
 101  C(K)=(C(K)+B(K)*C(K-1))*TERM
      DO 102 KK=KS,NT
      K=KS+NT-KK
 102  F(I,J,K)=A(K)*F(I,J,K+1)+C(K)
 100  CONTINUE
      RETURN
C-----BLOCK CORRECTION FOR PRESSURE FIELD
  2   CONTINUE
      DO 201 I=IS,LT
      A(I)=0.0
      B(I)=0.0
      C(I)=0.0
      D(I)=0.0
      DO 201 J=JS,MT
      DO 201 K=KS,NT
      A(I)=A(I)-AE(I,J,K)
      B(I)=B(I)-AW(I,J,K)
      C(I)=C(I)+SU(I,J,K)
 201  CONTINUE
      DO 202 I=IS,LT
      ANAB=-A(I)-B(I)
 202  D(I)=ANAB*1.5
      DO 203 I=IS+1,LT
      TERM=B(I)/D(I-1)
      D(I)=D(I)-TERM*A(I-1)
 203  C(I)=C(I)-TERM*C(I-1)
      F(LT,JS,KS)=C(LT)/D(LT)
      DO 204 II=IS+1,LT
      I=IS+LT-II
 204  F(I,JS,KS)=(C(I)-A(I)*F(I+1,JS,KS))/D(I)
      DO 205 I=IS,LT
      PPBLK=F(I,JS,KS)
      DO 205 J=JS+1,MT
      DO 205 K=KS+1,NT
      F(I,J,K)=PPBLK
 205  CONTINUE
      RETURN
      END
      SUBROUTINE BOUNC(IE,F)
      DIMENSION F(21,18,10)
      COMMON
     1/VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2  DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3  PP(21,18,10),W(21,18,10),TM(21,18,10)
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/PCCR/ DU(21,18,10),DV(21,18,10),DW(21,18,10)
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
     1/COEF/ AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2   SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3   AS(21,18,10),AT(21,18,10),AB(21,18,10),APO(21,18,10)
      COMMON
     1/TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJO(21,18,10),
     2   CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3   EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     3   SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
C-----EVALUATE WALL BOUNDARY CONDITIONS USING WALL FUNCTIONS
      DO 150 III=1,IITO
      I=IBC(III)
      J=JBC(III)
      K=KBC(III)
      GO TO (1,2,3,4,5,6), IITY(III)
  1   CONTINUE
C-----NORTH
      YP=YN(III)
      YP1=YN1(III)
      CDK=DK(I,J,K)
      CDE=DE(I,J,K)
      CDEN=DEN(I,J,K)
      CALL WALLFN(IE,YP,YP1,CDK,CDE,CDEN,SINX(III),SINY(III),
     1     SINZ(III),F(I,J+1,K),AN(I,J,K),SU(I,J,K),SP(I,J,K),
     2     TAUN(III),YPLN(III),GEN(I,J,K),VISE(I,J,K),
     3     U(I,J+1,K),V(I,J+1,K),W(I,J+1,K),U(I,J,K),V(I,J,K),
     4     W(I,J,K),U(I,J-1,K),V(I,J-1,K),W(I,J-1,K),
     5    -SPK(I,J,K),SPK(I,J,K),ANV1(III),ANW1(III),TJO(I,J,K))
      GO TO 150
  2   CONTINUE
C-----SOUTH
      YP=YN(III)
      YP1=YN1(III)
      CDK=DK(I,J,K)
      CDE=DE(I,J,K)
      CDEN=DEN(I,J,K)
      CALL WALLFN(IE,YP,YP1,CDK,CDE,CDEN,SINX(III),SINY(III),
     1     SINZ(III),F(I,J-1,K),AS(I,J,K),SU(I,J,K),SP(I,J,K),
     2     TAUN(III),YPLN(III),GEN(I,J,K),VISE(I,J,K),
     3     U(I,J-1,K),V(I,J-1,K),W(I,J-1,K),U(I,J,K),V(I,J,K),
     4     W(I,J,K),U(I,J+1,K),V(I,J+1,K),W(I,J+1,K),
     5    -SPK(I,J,K),SPK(I,J,K),ANV1(III),ANW1(III),TJO(I,J,K))
      GO TO 150
  3   CONTINUE
C-----EAST
      YP=YN(III)
      YP1=YN1(III)
      CDK=DK(I,J,K)
      CDE=DE(I,J,K)
      CDEN=DEN(I,J,K)
      CALL WALLFN(IE,YP,YP1,CDK,CDE,CDEN,SINX(III),SINY(III),
     1     SINZ(III),F(I+1,J,K),AE(I,J,K),SU(I,J,K),SP(I,J,K),
     2     TAUN(III),YPLN(III),GEN(I,J,K),VISE(I,J,K),
     3     U(I+1,J,K),V(I+1,J,K),W(I+1,J,K),U(I,J,K),V(I,J,K),
     4     W(I,J,K),U(I-1,J,K),V(I-1,J,K),W(I-1,J,K),
     5    -SPK(I,J,K),SPK(I,J,K),ANV1(III),ANW1(III),TJO(I,J,K))
      GO TO 150
  4   CONTINUE
C-----WEST
      YP=YN(III)
      YP1=YN1(III)
      CDK=DK(I,J,K)
      CDE=DE(I,J,K)
      CDEN=DEN(I,J,K)
      CALL WALLFN(IE,YP,YP1,CDK,CDE,CDEN,SINX(III),SINY(III),
     1     SINZ(III),F(I-1,J,K),AW(I,J,K),SU(I,J,K),SP(I,J,K),
     2     TAUN(III),YPLN(III),GEN(I,J,K),VISE(I,J,K),
     3     U(I-1,J,K),V(I-1,J,K),W(I-1,J,K),U(I,J,K),V(I,J,K),
     4     W(I,J,K),U(I+1,J,K),V(I+1,J,K),W(I+1,J,K),
     5    -SPK(I,J,K),SPK(I,J,K),ANV1(III),ANW1(III),TJO(I,J,K))
      GO TO 150
  5   CONTINUE
C-----TOP
      YP=YN(III)
      YP1=YN1(III)
      CDK=DK(I,J,K)
      CDE=DE(I,J,K)
      CDEN=DEN(I,J,K)
      CALL WALLFN(IE,YP,YP1,CDK,CDE,CDEN,SINX(III),SINY(III),
     1     SINZ(III),F(I,J,K+1),AT(I,J,K),SU(I,J,K),SP(I,J,K),
     2     TAUN(III),YPLN(III),GEN(I,J,K),VISE(I,J,K),
     3     U(I,J,K+1),V(I,J,K+1),W(I,J,K+1),U(I,J,K),V(I,J,K),
     4     W(I,J,K),U(I,J,K-1),V(I,J,K-1),W(I,J,K-1),
     5    -SPK(I,J,K),SPK(I,J,K),ANV1(III),ANW1(III),TJO(I,J,K))
      GO TO 150
  6   CONTINUE
C-----BOTTOM
      YP=YN(III)
      YP1=YN1(III)
      CDK=DK(I,J,K)
      CDE=DE(I,J,K)
      CDEN=DEN(I,J,K)
      CALL WALLFN(IE,YP,YP1,CDK,CDE,CDEN,SINX(III),SINY(III),
     1     SINZ(III),F(I,J,K-1),AB(I,J,K),SU(I,J,K),SP(I,J,K),
     2     TAUN(III),YPLN(III),GEN(I,J,K),VISE(I,J,K),
     3     U(I,J,K-1),V(I,J,K-1),W(I,J,K-1),U(I,J,K),V(I,J,K),
     4     W(I,J,K),U(I,J,K+1),V(I,J,K+1),W(I,J,K+1),
     5    -SPK(I,J,K),SPK(I,J,K),ANV1(III),ANW1(III),TJO(I,J,K))
 150  CONTINUE
      RETURN
      END
      SUBROUTINE WALLFN(IE,YP,YP1,DK1,DE1,DEN1,SINX1,SINY1,
     1     SINZ1,F1,AN1,SU1,SP1,
     2     TAUN1,YPLN1,GEN1,VISE1,
     3     U1,V1,W1,UO,VO,
     4     WO,U2,V2,W2,
     5     SUK1,SPK1,ANV2,ANW2,TJO1)
      COMMON
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
C-----WALL FUNCTIONS USING LOGARITHMIC WALL LAW
      GO TO (1,2,3,4,5,6,7), IE
  1   CONTINUE
C-----U
      SQRTK=SQRT(DK1)
      YPLN1=DEN1*SQRTK*CMU1*YP/VISC
      IF(YPLN1 .LE. 11.63) GO TO 111
      TMULT=DEN1*CMU1*SQRTK*CK/ALOG(E*YPLN1)
      GO TO 112
 111  TMULT=VISC/YP
 112  TAUN1=-TMULT
      PTA=TJO1*SINX1*TMULT/YP1
      SP1=SP1-PTA
      SU1=SU1+PTA*F1
      P1=SQRT(1.-SINX1**2)
      ANV2=AN1*SQRT(1.-SINY1**2)
      ANW2=AN1*SQRT(1.-SINZ1**2)
      AN1=P1*AN1
      RETURN
  2   CONTINUE
C-----V
      TMULT=-TAUN1
      PTA=TJO1*SINY1*TMULT/YP1
      SP1=SP1-PTA
      SU1=SU1+PTA*F1
      AN1=ANV2
      RETURN
  3   CONTINUE
C-----W
      TMULT=-TAUN1
      PTA=TJO1*SINZ1*TMULT/YP1
      SP1=SP1-PTA
      SU1=SU1+PTA*F1
      AN1=ANW2
      RETURN
  4   CONTINUE
C-----TM
      RETURN
  5   CONTINUE
C-----K
      SQRTK=SQRT(DK1)
      IF(YPLN1 .LE. 11.63) GO TO 511
      DITM=DEN1*CMU2*SQRTK*ALOG(E*YPLN1)/(CK*YP)
      GO TO 512
 511  DITM=DEN1*CMU2*SQRTK*YPLN1/YP
 512  CONTINUE
      DDU=UO-U1
      DDV=VO-V1
      DDW=WO-W1
      P1=DDU**2+DDV**2+DDW**2
      DDU=U2-U1
      DDV=V2-V1
      DDW=W2-W1
      P2=DDU**2+DDV**2+DDW**2
      YP2=2*YP1-YP
      GENR=GEN1-VISE1*P2/YP2/YP2
      IF(GENR .LE. 0.0) GENR=0.0
      SU1=SU1+TJO1*(P1*TAUN1**2/VISE1+GENR-GEN1)
      SP1=TJO1*(SPK1-DITM)
      AN1=0.0
      RETURN
  6   CONTINUE
C-----E
      TERM=CMU2/(CK*YP)
      SU1=HINUM*TERM*DK1**1.5
      SP1=-HINUM
      RETURN
  7   CONTINUE
C-----PF
      RETURN
      END
      SUBROUTINE SYMOUT(IC,IE,IS,IT,JS,JT,KS,KT)
      COMMON
     1/VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2  DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3  PP(21,18,10),W(21,18,10),TM(21,18,10)
     1/PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN
     1/PCCR/ DU(21,18,10),DV(21,18,10),DW(21,18,10)
     1/COEF/ AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2   SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3   AS(21,18,10),AT(21,18,10),AB(21,18,10),APO(21,18,10)
      COMMON
     1/TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJO(21,18,10),
     2   CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3   EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     3   SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1/UNSTDY/UO(21,18,10),VO(21,18,10),WO(21,18,10),DKO(21,18,10),
     2   DEO(21,18,10),DENO(21,18,10),TMO(21,18,10)
     1/LIMT/ L,M,LT,MT,L1,L2,M1,M2,LO,MO,ISWU,ISWV,ISWP,ISWK,ISWE,
     2   ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,NO,ISWW,IG,NT,ALC,DTT
C-----SYMMETRIC, CYCLIC AND EXIT CONDITIONS AND LINK MODIFICATIONS
      GO TO (1,2,3), IC
  1   CONTINUE
C-----BOTTOM
      K=1
      DO 10 I=1,L
      DO 10 J=2,MT
      U(I,J,K)=U(I,J,K+1)
      V(I,J,K)=V(I,J,K+1)
      W(I,J,K)=0.0
      TM(I,J,K)=TM(I,J,K+1)
      DK(I,J,K)=DK(I,J,K+1)
      DE(I,J,K)=DE(I,J,K+1)
  10  CONTINUE
C-----EAST OUT (BASED ON INFLOW MASS FLOW RATE)
      I=IT
      FLOW=0.0
      ARDEN=0.0
      DO 50 J=2,JT+1
      DO 50 K=2,KT+1
      UC=(V(I,J,K)+V(I,J-1,K)+V(I,J,K-1)+V(I,J-1,K-1))*0.25
      DENC=(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I,J,K-1)+DEN(I,J-1,K-1))*0.25
      P1=(X(I,J,K)+X(I,J,K-1)-X(I,J-1,K)-X(I,J-1,K-1))*0.5
      P2=(Y(I,J,K)+Y(I,J,K-1)-Y(I,J-1,K)-Y(I,J-1,K-1))*0.5
      P3=(Z(I,J,K)+Z(I,J,K-1)-Z(I,J-1,K)-Z(I,J-1,K-1))*0.5
      Q1=(X(I,J,K)+X(I,J-1,K)-X(I,J,K-1)-X(I,J-1,K-1))*0.5
      Q2=(Y(I,J,K)+Y(I,J-1,K)-Y(I,J,K-1)-Y(I,J-1,K-1))*0.5
      Q3=(Z(I,J,K)+Z(I,J-1,K)-Z(I,J,K-1)-Z(I,J-1,K-1))*0.5
      AREA=SQRT(P1*P1+P2*P2+P3*P3)*SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
      FLOW=FLOW+DENC*AREA*UC
      ARDEN=ARDEN+DENC*AREA
  50  CONTINUE
      UINC=(FLOW-FLOWIN)/ARDEN
      DO 60 J=2,JT
      DO 60 K=2,KT
      U(I+1,J,K)=U(I,J,K)
      V(I+1,J,K)=V(I,J,K)-UINC
      W(I+1,J,K)=W(I,J,K)
      TM(I+1,J,K)=TM(I,J,K)
      DK(I+1,J,K)=DK(I,J,K)
      DE(I+1,J,K)=DE(I,J,K)
  60  CONTINUE
      RETURN
C-----LINK COEFF. MODIFICATIONS
  2   CONTINUE
C-----EAST OUT
      I=IT
      DO 200 J=2,JT
      DO 200 K=2,KT
      AE(I,J,K)=0.0
 200  CONTINUE
      RETURN
C-----UPDATE UNSTEADY COEFF.
  3   CONTINUE
      IF(DTT .NE. 0.0) GO TO 301
      DO 300 I=IS,IT
      DO 300 J=JS,JT
      DO 300 K=KS,KT
 300  APO(I,J,K)=0.0
      RETURN
 301  CONTINUE
      DO 310 I=IS,IT
      DO 310 J=JS,JT
      DO 310 K=KS,KT
      APO(I,J,K)=DENO(I,J,K)/DTT
      UO(I,J,K)=U(I,J,K)
      VO(I,J,K)=V(I,J,K)
      WO(I,J,K)=W(I,J,K)
      TMO(I,J,K)=TM(I,J,K)
      DKO(I,J,K)=DK(I,J,K)
      DEO(I,J,K)=DE(I,J,K)
 310  CONTINUE
      RETURN
      END
      SUBROUTINE WALVAL(PW,IS,IT,JS,JT,KS,KT,F)
      DIMENSION F(21,18,10)
      COMMON
     1/TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2 YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3 YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4 TALW(800),GEN(21,18,10),MC(21,18,10),IJLO(21,18,10),IITO
C-----ASSIGN WALL VALUES
      DO 10 J=JS,JT
      DO 10 K=KS,KT
      F(IS-1,J,K)=PW*F(IS,J,K)
  10  F(IT+1,J,K)=PW*F(IT,J,K)
      DO 20 I=IS-1,IT+1
      DO 20 K=KS,KT
      F(I,JS-1,K)=PW*F(I,JS,K)
  20  F(I,JT+1,K)=PW*F(I,JT,K)
      DO 30 I=IS-1,IT+1
      DO 30 J=JS-1,JT+1
      F(I,J,KS-1)=PW*F(I,J,KS)
  30  F(I,J,KT+1)=PW*F(I,J,KT)
      DO 40 III=1,IITO
      I=IBC(III)
      J=JBC(III)
      K=KBC(III)
      GO TO (1,2,3,4,5,6), IITY(III)
  1   F(I,J+1,K)=PW*F(I,J,K)
      GO TO 40
  2   F(I,J-1,K)=PW*F(I,J,K)
      GO TO 40
  3   F(I+1,J,K)=PW*F(I,J,K)
      GO TO 40
  4   F(I-1,J,K)=PW*F(I,J,K)
      GO TO 40
  5   F(I,J,K+1)=PW*F(I,J,K)
      GO TO 40
  6   F(I,J,K-1)=PW*F(I,J,K)
  40  CONTINUE
      RETURN
      END
