      PROGRAM SKETCH
      CHARACTER*80 FN,OUT,T
      REAL X(5000),Y(5000)
      INTEGER IX(5000),IY(5000)
      CALL GETARG(1,FN)
      CALL GETARG(2,OUT)
      CALL GETARG(3,T)
      READ(T,*) FH,DY,WF
      U=0.84667
      OPEN(1,FILE=FN)
      READ(1,*) NW,NH,NL,NP,NT
      CALL METAFL('PDF')
      CALL SETFIL(OUT(1:LEN_TRIM(OUT)))
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(NW*U),NINT(NH*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      DO 20 L=1,NL
      READ(1,*) N,W
      READ(1,*) (X(I),Y(I),I=1,N)
      CALL LINWID(MAX(1,NINT(W*WF*U)))
      CALL STRTPT(X(1)*U,Y(1)*U)
      DO 10 I=2,N
  10  CALL CONNPT(X(I)*U,Y(I)*U)
  20  CONTINUE
      CALL LINWID(1)
      CALL SHDPAT(16)
      DO 40 L=1,NP
      READ(1,*) N
      READ(1,*) (X(I),Y(I),I=1,N)
      DO 30 I=1,N
      IX(I)=NINT(X(I)*U)
  30  IY(I)=NINT(Y(I)*U)
  40  CALL AREAF(IX,IY,N)
      CALL SHDPAT(0)
      CALL PSFONT('Helvetica-Bold')
      CALL TXTJUS('LEFT')
      CALL TXTJUS('BOTTOM')
      DO 60 L=1,NT
      READ(1,*) X0,X1,YB,CH
      READ(1,'(A)') T
      K=1
      IF(T(1:1).EQ.'@') THEN
        CALL PSFONT('Symbol')
        K=2
      ENDIF
      CALL HEIGHT(MAX(2,NINT(CH*FH*U)))
      CALL CHAWTH(1.)
      NL0=NLMESS(T(K:LEN_TRIM(T)))
      IF(NL0.GT.0) CALL CHAWTH(((X1-X0)*U)/NL0)
      CALL MESSAG(T(K:LEN_TRIM(T)),NINT(X0*U),NINT((YB+DY)*U))
      CALL PSFONT('Helvetica-Bold')
  60  CONTINUE
      CALL DISFIN
      END
