      PROGRAM FIG11
      REAL ZL(20),P(14,2)
      CHARACTER*1 PN(2)
      CHARACTER*12 FN
      DATA PN/'a','b'/
      U=0.84667
      OPEN(2,FILE='fig11_style.dat')
      READ(2,*) ((P(I,K),I=1,14),K=1,2),WL,WW,HL,CWL,ZM1,ZM2,ZM3
      CLOSE(2)
      CALL METAFL('PDF')
      CALL SETFIL('fig11_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(1442*U),NINT(990*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      CALL PSFONT('Helvetica-Bold')
      DO 50 K=1,2
      X0=P(1,K)
      Y0=P(2,K)
      SX=P(3,K)
      SY=P(4,K)
      ZL(1)=ZM1
      ZL(2)=ZM2
      DO 5 L=1,9
   5  ZL(2+L)=0.1*L
      ZL(12)=0.
      DO 6 L=1,4
   6  ZL(12+L)=P(4+L,K)
      ZL(17)=P(9,K)
      ZL(18)=ZM3
      CALL LINWID(MAX(1,NINT(WL*U)))
      WRITE(FN,'(A,A,A)') 'fig11',PN(K),'.dat'
      CALL CPANEL(FN,18,ZL,X0,Y0,SX,SY,0.,23.,0.,2.,U)
      CALL LINWID(MAX(1,NINT(WW*U)))
      YT=Y0-2*SY
      YS=Y0-SY
      XE=X0+23*SX
      CALL PL(X0,Y0,XE,Y0,U)
      CALL PL(X0,YT,XE,YT,U)
      CALL PL(X0-P(10,K),Y0,X0,Y0,U)
      CALL PL(X0-P(10,K),YS,X0,YS,U)
      CALL PL(X0-P(10,K),Y0,X0-P(10,K),YS,U)
      CALL PL(X0,Y0,X0,YS,U)
      CALL PL(XE,Y0+5,XE,Y0+20,U)
      YX=Y0+P(11,K)
      CALL PL(X0,YX-22,X0,YX+22,U)
      CALL PL(X0,YX,X0+60,YX,U)
      CALL HEAD(X0+75,YX,X0+60,YX,U)
      CALL CHAWTH(CWL)
      CALL HEIGHT(NINT(HL*U))
      CALL TXTJUS('LEFT')
      CALL TXTJUS('MIDDLE')
      CALL MESSAG('X',NINT((X0+92)*U),NINT(YX*U))
      CALL TXTJUS('CENT')
      CALL MESSAG('23h',NINT(XE*U),NINT((Y0+P(12,K))*U))
      CALL MESSAG('('//PN(K)//')',NINT(P(13,K)*U),NINT(P(14,K)*U))
  50  CONTINUE
      CALL DISFIN
      END
      SUBROUTINE PL(X1,Y1,X2,Y2,U)
      CALL LINE(NINT(X1*U),NINT(Y1*U),NINT(X2*U),NINT(Y2*U))
      RETURN
      END
      SUBROUTINE HEAD(XT,YT,XB,YB,U)
      INTEGER NX(3),NY(3)
      D=SQRT((XT-XB)**2+(YT-YB)**2)
      EX=(XT-XB)/D
      EY=(YT-YB)/D
      NX(1)=NINT(XT*U)
      NY(1)=NINT(YT*U)
      NX(2)=NINT((XB-EY*5)*U)
      NY(2)=NINT((YB+EX*5)*U)
      NX(3)=NINT((XB+EY*5)*U)
      NY(3)=NINT((YB-EX*5)*U)
      CALL SHDPAT(16)
      CALL AREAF(NX,NY,3)
      CALL SHDPAT(0)
      RETURN
      END
