      SUBROUTINE MESH5(TAG)
      CHARACTER*1 TAG
      REAL X(31,31),Y(31,31),XP(31),YP(31)
      U=0.84667
      OPEN(2,FILE='fig5'//TAG//'_style.dat')
      READ(2,*) BX0,BY0,BX1,BY1,W1,W2,RD
      CLOSE(2)
      OPEN(1,FILE='fig5'//TAG//'_nodes.dat')
      READ(1,*) NI,NJ
      DO 5 K=1,NI*NJ
      READ(1,*) I,J,A,B
      X(I,J)=BX0+A*(BX1-BX0)
   5  Y(I,J)=BY1-B*(BY1-BY0)
      CLOSE(1)
      CALL METAFL('PDF')
      CALL SETFIL('fig5'//TAG//'_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(1169*U),NINT(1146*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      CALL LINWID(MAX(1,NINT(W1*U)))
      DO 20 J=1,NJ
      DO 10 I=1,NI
      XP(I)=X(I,J)*U
  10  YP(I)=Y(I,J)*U
      CALL STRTPT(XP(1),YP(1))
      DO 15 I=2,NI
  15  CALL CONNPT(XP(I),YP(I))
  20  CONTINUE
      CALL LINWID(MAX(1,NINT(W2*U)))
      DO 40 I=1,NI
      CALL STRTPT(X(I,1)*U,Y(I,1)*U)
      DO 35 J=2,NJ
  35  CALL CONNPT(X(I,J)*U,Y(I,J)*U)
  40  CONTINUE
      IF(RD.GT.0.3) THEN
        CALL SHDPAT(16)
        DO 50 J=1,NJ
        DO 50 I=1,NI
  50    CALL CIRCLE(NINT(X(I,J)*U),NINT(Y(I,J)*U),MAX(1,NINT(RD*U)))
        CALL SHDPAT(0)
      ENDIF
      CALL DISFIN
      RETURN
      END
