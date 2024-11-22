      PROGRAM CNS3D
C     BY: Y. S. CHEN
C
C*****[ CURVILINEAR N-S CODE FOR 3-D INCOMPRESSIBLE FLOW ]**********
C
      COMMON
     1   /VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2     DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3     PP(21,18,10),W(21,18,10),TM(21,18,10),
     1   /PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN,
     1   /PCCR/ CU(21,18,10),DV(21,18,10),DW(21,18,10)
     1   /TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2     YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3     YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4     TALW(800),GEN(21,18,10),MC(21,18,10),IJL0(21,18,10),IIT0
     1   /COEF/AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2     SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3     AS(21,18,10),AT(21,18,10),AB(21,18,10),AP0(21,18,10)
      COMMON
     1     /TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJ0(21,18,10),
     2     CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3     EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     4     SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1     /LIMIT/ L,M,LT,MT,L1,L2,M1,M2,L0,M0ISWU,ISWV,ISWP,ISWK,ISWE,
     2     ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,N0,ISWW,IG,NT,ALC,DTT
      COMMON
     1     /TTRAN/TXXE(21,18,10),TXXW(21,18,10),TYYN(21,18,10),
     2     TYYS(21,18,10),TZZT(21,18,10),TZZB(21,18,10),
     3     TYXE(21,18,10),TYXW(21,18,10),TYZT(21,18,10),
     4     TYZS(21,18,10),TYXN(21,18,10),TXYS(21,18,10),
     5     TXZT(21,18,10),TXZB(21,18,10),TZXE(21,18,10),
     6     TZXW(21,18,10),TZYN(21,18,10),TZYS(21,18,10)
      COMMON
     1     /UNSTDY/U0(21,18,10),V0(21,18,10),W0(21,18,10),DK0(21,18,10),
     2     DE0(21,18,10),DEN0(21,18,10),TM0(21,18,10)
      LOGICAL INS0U,INS0V,INS0P,INS0K,INS0E,INSR0,INS0W,INS0T
C*****[ INPUT DATA GUIDE ]*********************************************
C  NLIMT : MAXIMUM NO. OF ITERATIONS LIMIT
C
C  IG    = 1: LAMINAR
C          2: TURBULENT (K-E MODEL)
C
C  ISWP  : NO. OF SWEEPS FOR SOLVING THE P' EQUATION (PP).
C
C**********************************************************************
C-----INPUT DATA (PROBLEM CONTROL SETTING)
      READ(5,100) NLIMT,IG,ISWP,ITT
      READ(5,200) ALU,ALV,ALW,ALP,ALK,ALE,ALVIS,ALC
      READ(5,200) RENL,DTT
C-----CONSTANTS
      EREXT=1.0E-3
      ISWL=7
      ISWV=7
      ISWW=7
      ISWK=5
      ISWE=5
      DENIN=1.0
      VISC=1.0/RENL
      SIGL=1.0
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
      INS0U=.TRUE.
      INS0V=.TRUE.
      INS0P=.TRUE.
      CMU1=CMU**0.25
      CMU2=CMU**0.75
      INS0U=.TRUE.
      INS0V=.TRUE.
      INS0P=.TRUE.
      INS0K=.TRUE.
      INS0E=.TRUE.
      INS0P=.TRUE.
      INS0W=.TRUE.
      INS0T=.TRUE.
      IF(IG .EQ. 2) GO TO 10
      INS0K=.TRUE.
      INS0E=.TRUE.
      INSR0=.TRUE.
 10   CONTINUE
C*****[READ IN INITIAL FLOW FIELDS FROM RESTART FILE (LU = 8)]********
      READ(8,100) L,M,N,L1,L2,M1,M2,N1,N2
      L0=L+1
      M0=M+1
      N0=N+1
      LT=L-1
      MT=M-1
      NT=N-1
C-----INITILIZE VARIABLES
      CALL INIT
C-----RESTART FILE
      DO 50 K=1,N
         DO 50 I=1,L
            READ(8,400)
            DO 50 J=1,M
               READ(8,500) X(I,J,K),Y(I,J,K),Z(I,J,K),U(I,J,K),V(I,J,K),
     1              W(I,J,K),P(I,J,K),TM(I,J,K),DK(I,J,K),DE(I,J,K),
     2              VISCTM,DENTM
               U0(I,J,K)=U(I,J,K)
               V0(I,J,K)=V(I,J,K)
               TM0(I,J,K)=TM(I,J,K)
               DE0(I,J,K)=DE(I,J,K)
               DEN0(I,J,K)=DEN(I,J,K)
 50         CONTINUE
C************************************************************
C--------GET BOUNDARY CONTROL PARAMETRS
      CALL DIRCOS
C--------SET BOUNDARY TURBULENCE PARAMETRS TO ZERO
            DO 121 I=1,L
               DO 121 J=1,M
                  DO 121 K=1,N
                     IF (MC(I,J,K) .NE. 0) GO TO 122
                     GO TO 121
 122                 DK(I,J,K) = 0.0
                     DE(I,J,K) = 0.0
                     U(I,J,K) = 0.0
                     V(I,J,K) = 0.0
                     W(I,J,K) = 0.0
 121              CONTINUE
C--------CALCULATE GRID TRANFORAMTION COEFFICIENTS
      CALL TRANF
C--------TURBULENT VISCOSITY
      IF(INSR0) CALL NEWVIS
C--------CALCULATE INLET MASS FLOW RATE
      FLOWIN=0.0
      I=1
      DO 45 J=2,M
         DO 45 K=2,N
            UC=(U(I,J,K)+U(I,J-1,K)+U(I,J,K-1)+U(I,J-1,K-1))*0.25
            DENC=(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I,J,K-1)+DEN(I,J-1,K-1))
     1           * 0.25
            P1=(X(I,J,K)+X(I,J,K-1)-X(I,J-1,K)-X(I,J-1,K-1))*0.5
            P2=(Y(I,J,K)+Y(I,J,K-1)-Y(I,J-1,K)-Y(I,J-1,K-1))*0.5
            P3=(Z(I,J,K)+Z(I,J,K-1)-Z(I,J-1,K)-Z(I,J-1,K-1))*0.5
            Q1=(X(I,J,K)+X(I,J-1,K)-X(I,J,K-1)-X(I,J-1,K-1))*0.5
            Q2=(Y(I,J,K)+Y(I,J-1,K)-Y(I,J,K-1)-Y(I,J-1,K-1))*0.5
            Q3=(Z(I,J,K)+Z(I,J-1,K)-Z(I,J,K-1)-Z(I,J-1,K-1))*0.5
            AREA=SQRT(P1*P1+P2*P2+P3*P3)*SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
            FLOWIN=FLOWIN+DENC*AREA*UC
 45      CONTINUE
      ITO=1
C--------TRANSIENT PROCESS START
 2    CONTINUE
      CALL SYMOUT(3,1,2,L,2,M,2,N)
      ITER=1
C--------SOLUTION PROCEDURES START
 1    CONTINUE
      CALL SYMOUT(1,1,2,LT,2,MT,2,NT)
      IF(INS0U) CALL SOLVEQ(1,ISWU,ALU,SIGU,ERRU,U,U0)
      IF(INS0V) CALL SOLVEQ(2,ISWV,ALV,SIGU,ERRV,V,V0)
      IF(INS0W) CALL SOLVEQ(3,ISWW,ALW,SIGU,ERRW,W,W0)
      IF(INS0T) CALL SOLVEQ(4,ISWW,ALW,SIGU,ERRW,TM,TM0)
      IF(INS0K) CALL SOLVEQ(5,ISWK,ALK,SIGK,ERRK,DK,DK0)
      IF(INS0E) CALL SOLVEQ(6,ISWE,ALE,SIGE,ERRE,DE,DE0)
      IF(INS0P) CALL SOLVEQ(0,ISWP,ALP,SIGU,ERRM,PP,PP)
      IF(INSR0) CALL NEWVIS
C--------CONVERGENCE CHECK
      WRITE(6,300) ITER,ERRU,ERRV,ERRW,ERRM,ERRK,ERRE,U(7,2,6)
      ERRMAX=ERRM+ERRU+ERRV+ERRW
      IF(ITER .GE. 20 .AND. ERRMAX .GT. 1.E03) GO TO 99
      IF(ITER .GE. NLIMT .OR. ERRMAX .LE. EREXT) GO TO 99
      ITER=ITER+1
      GO TO 1
C--------PRINT OUT SOLUTION
 99   CONTINUE
      WRITE(7,100) L,M,N,L1,L2,M1,M3,N1,N2
      DO 901 K=1,N
         DO 901 I=1,L
            WRITE(7,400)
            DO 902 J=1,M
               XV=X(I,J,K)
               YV=Y(I,J,K)
               ZV=Z(I,J,K)
               WRITE(7,500)XV,YV,ZV,U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),
     1              TM(I,J,K),DK(I,J,K),DE(I,J,K),VISE(I,J,K),DEN(I,J,K)
 902        CONTINUE
 901     CONTINUE
      TIMT=DTT*IT0
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
     1   /TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2     YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3     YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4     TALW(800),GEN(21,18,10),MC(21,18,10),IJL0(21,18,10),IIT0
     1     /TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJ0(21,18,10),
     2     CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3     EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     4     SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1     /LIMIT/ L,M,LT,MT,L1,L2,M1,M2,L0,M0ISWU,ISWV,ISWP,ISWK,ISWE,
     2     ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,N0,ISWW,IG,NT,ALC,DTT
C-----SET DOMAIN BLOCKAGE CONTROL PARAMETERS
C-----SCALAR BLOCKAGE :   MC(I,J,K)=1
C-----PRESSURE BLOCKAGE : MC(I,J,K)=1
      DO 10 I=1,L0
         DO 10 J=1,M0
            DO 10 K=1,N0
               IF(J.EQ.1.OR.J.EQ.M.OR.K.EQ.N) MC(I,J,K)=1
C     TODO
C-----ADD BLOCKAGE AS NEEDED HERE
 10   CONTINUE
C-----CALCULATE BOUNDARY GRID SIZES AND ORIENTATIONS
      III=1
      DO 30 I=2,LT
         DO 30 J=2,MT
            DO 30 K=2,NT
 30   CONTINUE
      IIT0=III-1
      WRITE(6,100) L0,M0,N0,IIT0
 100  FORMAT(4I5)
      RETURN
      END

      SUBROUTINE TRANF
      RETURN
      END

      SUBROUTINE INIT
      RETURN
      END

      SUBROUTINE NEWVIS
      RETURN
      END

      SUBROUTINE SOLVEQ(IE,ISWF,ALF,SIGF,ERF,F,F0)
      DIMENSION F(21,18,10),F1(21,18,10),F0(21,18,10)
      RETURN
      END

      SUBROUTINE LINERX
      RETURN
      END

      SUBROUTINE BOUND(IE,F)
      RETURN
      END

      SUBROUTINE WALLFN
      RETURN
      END

      SUBROUTINE SYMOUT(IC,IE,IS,IT,JS,JT,KS,KT)
      COMMON
     1   /VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2     DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3     PP(21,18,10),W(21,18,10),TM(21,18,10),
     1   /PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN,
     1   /PCCR/ CU(21,18,10),DV(21,18,10),DW(21,18,10)
     1   /COEF/AP(21,18,10),SU(21,18,10),SP(21,18,10),SUK(21,18,10),
     2     SPK(21,18,10),AE(21,18,10),AW(21,18,10),AN(21,18,10),
     3     AS(21,18,10),AT(21,18,10),AB(21,18,10),AP0(21,18,10)
      COMMON
     1     /TRAN/ X(21,18,10),Y(21,18,10),Z(21,18,10),TJ0(21,18,10),
     2     CX(21,18,10),CY(21,18,10),CZ(21,18,10),
     3     EX(21,18,10),EY(21,18,10),EZ(21,18,10),
     4     SX(21,18,10),SY(21,18,10),SZ(21,18,10)
     1     /UNSTDY/U0(21,18,10),V0(21,18,10),W0(21,18,10),DK0(21,18,10),
     2     DE0(21,18,10),DEN0(21,18,10),TM0(21,18,10)
     1     /LIMIT/ L,M,LT,MT,L1,L2,M1,M2,L0,M0ISWU,ISWV,ISWP,ISWK,ISWE,
     2     ALU,ALV,ALP,ALK,ALE,ALVIS,ALW,N,N1,N2,N0,ISWW,IG,NT,ALC,DTT
C-----SYMMETRIC, CYCLIC AND EXIT CONDTIONS AND LINK MODIFICATIONS
      GO TO (1,2,3), IC
 1    CONTINUE
C-----BOTTOM
      K=1
      DO 10 I=1,L
         DO 10 J=2,MT
            U(I,J,K)=U(I,J,K+1)
            V(I,J,K)=V(I,J,K+1)
            W(I,J,K)=0
            TM(I,J,K)=TM(I,J,K+1)
            DK(I,J,K)=DK(I,J,K+1)
            DE(I,J,K)=DE(I,J,K+1)
 10   CONTINUE
C-----EAST OUT (BASED ON INFLOW MASS FLOW RATE)
      I=IT
      FLOW=0.0
      ARDEN=0.0
      DO 50 J=2,JT
         DO 50 K=2,KT
            UC=(V(I,J,K)+V(I,J-1,K)+V(I,J,K-1)+V(I,J-1,K-1))*0.25
            DENC=(DEN(I,J,K)+DEN(I,J-1,K)+DEN(I,J,K-1)+DEN(I,J-1,K-1))
     1           * 0.25
            P1=(X(I,J,K)+X(I,J,K-1)-X(I,J-1,K)-X(I,J-1,K-1))*0.5
            P2=(Y(I,J,K)+Y(I,J,K-1)-Y(I,J-1,K)-Y(I,J-1,K-1))*0.5
            P3=(Z(I,J,K)+Z(I,J,K-1)-Z(I,J-1,K)-Z(I,J-1,K-1))*0.5
            Q1=(X(I,J,K)+X(I,J-1,K)-X(I,J,K-1)-X(I,J-1,K-1))*0.5
            Q2=(Y(I,J,K)+Y(I,J-1,K)-Y(I,J,K-1)-Y(I,J-1,K-1))*0.5
            Q3=(Z(I,J,K)+Z(I,J-1,K)-Z(I,J,K-1)-Z(I,J-1,K-1))*0.5
            AREA=SQRT(P1*P1+P2*P2+P3*P3)*SQRT(Q1*Q1+Q2*Q2+Q3*Q3)
            FLOW=FLOW+DENC*AREA*UC
            ARDEN=ARDEN+DENC*AREA
 50   CONTINUE
      UINC=(FLOW-FLOWIN)/ARDEN
      DO 60 J=2,JT
         DO 60 K=2,KT
            U(I+1,J,K)=U(I,J,K)
            V(I+1,J,K)=V(I,J,K)-UINC
            W(I+1,J,K)=W(I,J,K)
            DK(I+1,J,K)=DK(I,J,K)
            DE(I+1,J,K)=DE(I,J,K)
 60   CONTINUE
      RETURN
C-----LINK COEF. MODIFICATIONS
 2    CONTINUE
C-----EAST OUT
      I=IT
      DO 200 J=2,JT
         DO 200 K=2,KT
            AE(I,J,K)=0.0
 200     CONTINUE
      RETURN
C-----UPDATE UNSTEADY COEFF.
 3    CONTINUE
      IF(DTT .NE. 0.0) GO TO 301
      DO 300 I=IS,IT
         DO 300 J=JS,JT
            DO 300 K=KS,KT
 300           AP0(I,J,K)=0.0
               RETURN
 301  CONTINUE
      DO 310 I=IS,IT
         DO 310 J=JS,JT
            DO 310 K=KS,KT
               AP0(I,J,K)=DEN0(I,J,K)/DTT
               U0(I,J,K)=U(I,J,K)
               V0(I,J,K)=V(I,J,K)
               W0(I,J,K)=W(I,J,K)
               TM0(I,J,K)=TM(I,J,K)
               DK0(I,J,K)=DK(I,J,K)
               DE0(I,J,K)=DE(I,J,K)
 310  CONTINUE
      RETURN
      END

      SUBROUTINE WALVAL(PW,IS,IT,JS,JTKS,KT,F)
      DIMENSION F(21,18,10)
      COMMON
     1   /TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2     YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3     YPLN(800),TAUN(800),IBC(800),JBC(800),KBC(800),IITY(800),
     4     TALW(800),GEN(21,18,10),MC(21,18,10),IJL0(21,18,10),IIT0
C-----ASSINGN WALL VALUES
      DO 10 J=JS,JT
         DO 10 K=KS,KT
            F(IS-1,J,K)=PW*F(IS,J,K)
 10         F(IT+1,J,K)=PW*F(IS,J,K)
      DO 20 I=IS-1,IT+1
         DO 20 K=KS,KT
            F(I,JS-1,K)=PW*F(I,JS,K)
 20         F(I,JT+1,K)=PW*F(I,JT,K)
      DO 30 I=IS-1,IT+1
         DO 30 J=JS-1,JT+1
            F(I,J,KS-1)=PW*F(I,J,KS)
 30         F(I,J,KT+1)=PW*F(I,J,KT)
      DO 40 III=1,IIT0
         I=IBC(III)
         J=JBC(III)
         K=KBC(III)
         GO TO (1,2,3,4,5,6), IITY(III)
 1       F(I,J+1,K)=PW*F(I,J,K)
         GO TO 40
 2       F(I,J-1,K)=PW*F(I,J,K)
         GO TO 40
 3       F(I+1,J,K)=PW*F(I,J,K)
         GO TO 40
 4       F(I-1,J,K)=PW*F(I,J,K)
         GO TO 40
 5       F(I,J,K+1)=PW*F(I,J,K)
         GO TO 40
 6       F(I,J,K-1)=PW*F(I,J,K)
 40   CONTINUE
      RETURN
      END
