      COMMON
     1   /VAR/U(21,18,10),V(21,18,10),P(21,18,10),DK(21,18,10),
     2     DE(21,18,10),ERRU,ERRV,ERRM,ERRK,ERRE,ERRW,
     3     PP(21,18,10),W(21,18,10),TM(21,18,10),
     1   /PRCP/ VISE(21,18,10),DEN(21,18,10),VISC,DENIN,FLOWIN,
     1   /PCCR/ CU(21,18,10),DV(21,18,10),DW(21,18,10)
     1   /TUR/ SIGK,SIGE,CMU,C1,C2,CMU1,CMU2,E,CK,HINUM,SMNUM,ANV1(800),
     2     YN(800),YN1(800),SINX(800),SINY(800),SINZ(800),ANW1(800),
     3     YPLN(800),TAUN(800),IBC(800),JBC(800),KBZ(800),IITY(800),
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
      LOGICAL INS0U,INS0V,INS0P,INS0K,INS0E,INPR0,INS0W,INS0T
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
      INPR0=.TRUE.
      INS0W=.TRUE.
      INS0T=.TRUE.
      IF(IG .EQ. 2) GO TO 10
      INS0K=.TRUE.
      INS0E=.TRUE.
      INPR0=.TRUE.
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
C      CALL INIT
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
C     CALL DIRCOS
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
 999  CONTINUE
 100  FORMAT(9I5)
 200  FORMAT(11F7.4)
 300  FORMAT(1X,I5,7E10.2)
 400  FORMAT(//)
 500  FORMAT(3F8.4,3E12.4,2X,6E11.4)
      STOP
      END
