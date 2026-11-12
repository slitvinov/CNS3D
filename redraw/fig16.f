      PROGRAM FIG16
      REAL X(40,20),Y(40,20),ZE(20),ET(20)
      U=0.84667
      OPEN(2,FILE='fig16_style.dat')
      READ(2,*) CX,CY,D,Q,XB0,YB0,XB1,YB1,W1,W2,HL,CWL,XL,YL
      CLOSE(2)
      OPEN(1,FILE='fig16_nodes.dat')
      READ(1,*) NI,NK
      DO 5 L=1,NI*NK
      READ(1,*) I,K,A,B
      X(I,K)=(CX+Q*D*A)*U
   5  Y(I,K)=(CY-D*B)*U
      CLOSE(1)
      OPEN(1,FILE='fig16_xsec.dat')
      READ(1,*) NZ,NR
      READ(1,*) (ZE(K),K=1,NZ)
      READ(1,*) (ET(K),K=1,NR)
      CLOSE(1)
      CALL METAFL('PDF')
      CALL SETFIL('fig16_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(2043*U),NINT(1713*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      CALL LINWID(MAX(1,NINT(W1*U)))
      DO 20 K=1,NK
      CALL STRTPT(X(1,K),Y(1,K))
      DO 15 I=2,NI
  15  CALL CONNPT(X(I,K),Y(I,K))
  20  CONTINUE
      CALL LINWID(MAX(1,NINT(W2*U)))
      DO 30 I=1,NI
  30  CALL PLINE(X(I,1),Y(I,1),X(I,NK),Y(I,NK))
      CALL LINWID(MAX(1,NINT(W1*U)))
      DO 40 K=1,NZ
      XX=(XB0+ZE(K)*(XB1-XB0))*U
  40  CALL PLINE(XX,YB0*U,XX,YB1*U)
      DO 45 K=1,NR
      YY=(YB0+ET(K)*(YB1-YB0))*U
  45  CALL PLINE(XB0*U,YY,XB1*U,YY)
      CALL PSFONT('Helvetica-Bold')
      CALL CHAWTH(CWL)
      CALL HEIGHT(NINT(HL*U))
      CALL TXTJUS('CENT')
      CALL TXTJUS('MIDDLE')
      CALL MESSAG('(b)',NINT(XL*U),NINT(YL*U))
      CALL DISFIN
      END
      SUBROUTINE PLINE(X1,Y1,X2,Y2)
      CALL LINE(NINT(X1),NINT(Y1),NINT(X2),NINT(Y2))
      RETURN
      END
