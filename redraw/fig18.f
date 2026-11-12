      PROGRAM FIG18
      REAL Z(200),R(200),WW(200),VR(200),B(4,4)
      CHARACTER*1 PN(4)
      CHARACTER*12 FN
      DATA PN/'a','b','c','d'/
      U=0.84667
      OPEN(2,FILE='fig18_style.dat')
      READ(2,*) ((B(I,K),I=1,4),K=1,4),SC,AL,AW,WV,WF,HL,CWL,DLB,AXL
      CLOSE(2)
      CALL METAFL('PDF')
      CALL SETFIL('fig18_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(2317*U),NINT(1740*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      CALL PSFONT('Helvetica-Bold')
      DO 50 K=1,4
      X0=B(1,K)
      Y0=B(2,K)
      X1=B(3,K)
      Y1=B(4,K)
      H=Y1-Y0
      CALL LINWID(MAX(1,NINT(WF*U)))
      CALL PLINE(X0*U,Y0*U,X1*U,Y0*U)
      CALL PLINE(X1*U,Y0*U,X1*U,Y1*U)
      CALL PLINE(X1*U,Y1*U,X0*U,Y1*U)
      CALL PLINE(X0*U,Y1*U,X0*U,Y0*U)
      WRITE(FN,'(A,A,A)') 'fig18',PN(K),'.dat'
      OPEN(1,FILE=FN)
      READ(1,*) N
      DO 10 I=1,N
  10  READ(1,*) Z(I),R(I),WW(I),VR(I)
      CLOSE(1)
      CALL LINWID(MAX(1,NINT(WV*U)))
      DO 20 I=1,N
      XA=X0+Z(I)/0.5*(X1-X0)
      YA=Y0+R(I)*H
  20  CALL ARROW(XA,YA,XA+SC*WW(I)*H,YA+SC*VR(I)*H,AL,AW,U)
      YZ=Y1+DLB+40
      CALL LINWID(MAX(1,NINT(WF*U)))
      CALL PLINE(X0*U,(YZ-28)*U,X0*U,(YZ+28)*U)
      CALL ARROW(X0,YZ,X0+AXL,YZ,AL*1.6,AW*1.6,U)
      CALL CHAWTH(CWL)
      CALL HEIGHT(NINT(HL*U))
      CALL TXTJUS('LEFT')
      CALL TXTJUS('MIDDLE')
      CALL MESSAG('z',NINT((X0+AXL+12)*U),NINT(YZ*U))
      CALL TXTJUS('CENT')
      CALL MESSAG('('//PN(K)//')',NINT((X0+X1)/2*U),NINT((Y1+DLB)*U))
  50  CONTINUE
      CALL DISFIN
      END
      SUBROUTINE PLINE(X1,Y1,X2,Y2)
      CALL LINE(NINT(X1),NINT(Y1),NINT(X2),NINT(Y2))
      RETURN
      END
      SUBROUTINE ARROW(X1,Y1,X2,Y2,AL,AW,U)
      INTEGER NX(3),NY(3)
      D=SQRT((X2-X1)**2+(Y2-Y1)**2)
      IF(D.LT.0.5) RETURN
      CALL LINE(NINT(X1*U),NINT(Y1*U),NINT(X2*U),NINT(Y2*U))
      IF(D.LT.AL) RETURN
      EX=(X2-X1)/D
      EY=(Y2-Y1)/D
      NX(1)=NINT(X2*U)
      NY(1)=NINT(Y2*U)
      NX(2)=NINT((X2-EX*AL-EY*AW)*U)
      NY(2)=NINT((Y2-EY*AL+EX*AW)*U)
      NX(3)=NINT((X2-EX*AL+EY*AW)*U)
      NY(3)=NINT((Y2-EY*AL-EX*AW)*U)
      CALL SHDPAT(16)
      CALL AREAF(NX,NY,3)
      CALL SHDPAT(0)
      RETURN
      END
