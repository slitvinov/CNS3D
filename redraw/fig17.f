      PROGRAM FIG17
      REAL X(500),Y(500),UU(500),VV(500),WX(100),WY(100),PP(4,3)
      REAL XL(3),YL(3)
      CHARACTER*1 PN(3)
      CHARACTER*3 LB
      CHARACTER*12 FN
      DATA PN/'a','b','c'/
      U=0.84667
      OPEN(2,FILE='fig17_style.dat')
      READ(2,*) ((PP(I,K),I=1,4),K=1,3),SC,AL,AW,WV,WF,HL,CWL,
     1  (XL(K),YL(K),K=1,3)
      CLOSE(2)
      OPEN(1,FILE='fig17_walls.dat')
      READ(1,*) NW
      DO 5 I=1,2*NW
   5  READ(1,*) WX(I),WY(I)
      CLOSE(1)
      CALL METAFL('PDF')
      CALL SETFIL('fig17_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(2548*U),NINT(1835*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      CALL PSFONT('Helvetica-Bold')
      DO 50 K=1,3
      CX=PP(1,K)
      CY=PP(2,K)
      D=PP(3,K)
      Q=PP(4,K)
      CALL LINWID(MAX(1,NINT(WF*U)))
      DO 10 IW=0,1
      CALL STRTPT((CX+Q*D*WX(IW*NW+1))*U,(CY-D*WY(IW*NW+1))*U)
      DO 10 I=2,NW
  10  CALL CONNPT((CX+Q*D*WX(IW*NW+I))*U,(CY-D*WY(IW*NW+I))*U)
      WRITE(FN,'(A,A,A)') 'fig17',PN(K),'.dat'
      OPEN(1,FILE=FN)
      READ(1,*) N
      DO 15 I=1,N
  15  READ(1,*) X(I),Y(I),UU(I),VV(I)
      CLOSE(1)
      CALL LINWID(MAX(1,NINT(WV*U)))
      DO 20 I=1,N
      XA=CX+Q*D*X(I)
      YA=CY-D*Y(I)
      XB=XA+SC*Q*D*UU(I)
      YB=YA-SC*D*VV(I)
  20  CALL ARROW(XA,YA,XB,YB,AL,AW,U)
      CALL CHAWTH(CWL)
      CALL HEIGHT(NINT(HL*U))
      CALL TXTJUS('CENT')
      CALL TXTJUS('MIDDLE')
      LB='('//PN(K)//')'
      CALL MESSAG(LB,NINT(XL(K)*U),NINT(YL(K)*U))
  50  CONTINUE
      CALL DISFIN
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
