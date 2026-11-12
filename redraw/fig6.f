      PROGRAM FIG6
      REAL X(1000),Y(1000),UU(1000),VV(1000),BX(4,2)
      CHARACTER*1 PN(2)
      CHARACTER*12 FN
      DATA PN/'a','b'/
      U=0.84667
      OPEN(2,FILE='fig6_style.dat')
      READ(2,*) ((BX(I,K),I=1,4),K=1,2),SC,RD,WV,WF,HL,CWL,XLA,XLB,YL
      CLOSE(2)
      CALL METAFL('PDF')
      CALL SETFIL('fig6_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(1910*U),NINT(1002*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      DO 50 K=1,2
      X0=BX(1,K)
      Y0=BX(2,K)
      X1=BX(3,K)
      Y1=BX(4,K)
      CALL LINWID(MAX(1,NINT(WF*U)))
      CALL PLINE(X0*U,Y0*U,X1*U,Y0*U)
      CALL PLINE(X1*U,Y0*U,X1*U,Y1*U)
      CALL PLINE(X1*U,Y1*U,X0*U,Y1*U)
      CALL PLINE(X0*U,Y1*U,X0*U,Y0*U)
      WRITE(FN,'(A,A,A)') 'fig6',PN(K),'.dat'
      OPEN(1,FILE=FN)
      READ(1,*) N
      DO 10 I=1,N
  10  READ(1,*) X(I),Y(I),UU(I),VV(I)
      CLOSE(1)
      CALL LINWID(MAX(1,NINT(WV*U)))
      CALL SHDPAT(16)
      DO 20 I=1,N
      IF(X(I).LE.0..OR.X(I).GE.1..OR.Y(I).LE.0..OR.Y(I).GE.1.) GOTO 20
      XA=X0+X(I)*(X1-X0)
      YA=Y1-Y(I)*(Y1-Y0)
      XB=XA+SC*UU(I)*(X1-X0)
      YB=YA-SC*VV(I)*(Y1-Y0)
      CALL PLINE(XA*U,YA*U,XB*U,YB*U)
      ALEN=SQRT((XB-XA)**2+(YB-YA)**2)
      IF(RD.GT.0.3.AND.ALEN.GT.3*RD)
     1  CALL CIRCLE(NINT(XB*U),NINT(YB*U),MAX(1,NINT(RD*U)))
  20  CONTINUE
      CALL SHDPAT(0)
  50  CONTINUE
      CALL PSFONT('Helvetica-Bold')
      CALL CHAWTH(CWL)
      CALL HEIGHT(NINT(HL*U))
      CALL TXTJUS('CENT')
      CALL TXTJUS('MIDDLE')
      CALL MESSAG('(a)',NINT(XLA*U),NINT(YL*U))
      CALL MESSAG('(b)',NINT(XLB*U),NINT(YL*U))
      CALL DISFIN
      END
      SUBROUTINE PLINE(X1,Y1,X2,Y2)
      CALL LINE(NINT(X1),NINT(Y1),NINT(X2),NINT(Y2))
      RETURN
      END
