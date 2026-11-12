      PROGRAM FIG15
      REAL ZL(6),XT(5),YT(5)
      INTEGER IP(2)
      CHARACTER*4 LB(5)
      DATA ZL/0.005,0.01,0.02,0.03,0.04,0.045/
      DATA LB/'.045','.04','.03','.02','.01'/
      U=0.84667
      OPEN(2,FILE='fig15_style.dat')
      READ(2,*) SXD,X0,Y0,SX,SY,WL,WW,HL,CWL,D1,G1,(XT(L),YT(L),L=1,5)
      CLOSE(2)
      CALL METAFL('PDF')
      CALL SETFIL('fig15_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(1877*U),NINT(435*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      CALL PSFONT('Helvetica-Bold')
      CALL LINWID(MAX(1,NINT(WL*U)))
      IP(1)=MAX(1,NINT(D1*U))
      IP(2)=MAX(1,NINT(G1*U))
      CALL MYLINE(IP,2)
      CALL CPANEL('fig15u.dat',6,ZL,X0,Y0,SX,SY,-4.,0.,0.,3.,U)
      XS=X0+4*SX
      XD=(1877.-XS)/SXD
      CALL CPANEL('fig15d.dat',6,ZL,XS,Y0,SXD,SY,0.,XD,0.,3.,U)
      CALL SOLID
      CALL LINWID(MAX(1,NINT(WW*U)))
      XS=X0+4*SX
      XE=XS+XD*SXD
      CALL PL(X0,Y0,XE,Y0,U)
      CALL PL(X0,Y0-3*SY,XE,Y0-3*SY,U)
      CALL PL(X0,Y0,X0,Y0-3*SY,U)
      CALL PL(X0,Y0-SY,XS,Y0-SY,U)
      CALL PL(XS,Y0,XS,Y0-SY,U)
      CALL CHAWTH(CWL)
      CALL HEIGHT(NINT(HL*U))
      CALL TXTJUS('CENT')
      CALL TXTJUS('MIDDLE')
      DO 10 L=1,5
  10  CALL MESSAG(LB(L),NINT(XT(L)*U),NINT(YT(L)*U))
      CALL DISFIN
      END
      SUBROUTINE PL(X1,Y1,X2,Y2,U)
      CALL LINE(NINT(X1*U),NINT(Y1*U),NINT(X2*U),NINT(Y2*U))
      RETURN
      END
