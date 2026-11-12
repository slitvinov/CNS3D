      PROGRAM FIG14
      REAL ZL(20)
      CHARACTER*4 LB
      U=0.84667
      OPEN(2,FILE='fig14_style.dat')
      READ(2,*) SXD,X0,Y0,SX,SY,ZM,ZN1,ZN2,ZN3,WL,WW,HL,CWL,XR,DLB
      CLOSE(2)
      CALL METAFL('PDF')
      CALL SETFIL('fig14_dislin.pdf')
      CALL FILMOD('DELETE')
      CALL PAGE(NINT(1879*U),NINT(492*U))
      CALL DISINI
      CALL ERRMOD('ALL','OFF')
      CALL PSFONT('Helvetica-Bold')
      DO 5 L=1,9
   5  ZL(L)=0.1*L
      ZL(10)=ZM
      ZL(11)=0.
      ZL(12)=ZN1
      ZL(13)=ZN2
      ZL(14)=ZN3
      CALL LINWID(MAX(1,NINT(WL*U)))
      CALL CPANEL('fig14u.dat',14,ZL,X0,Y0,SX,SY,-4.,0.,0.,3.,U)
      XS=X0+4*SX
      XD=(1879.-XS)/SXD
      CALL CPANEL('fig14d.dat',14,ZL,XS,Y0,SXD,SY,0.,XD,0.,3.,U)
      CALL LINWID(MAX(1,NINT(WW*U)))
      XS=X0+4*SX
      XE=XS+XD*SXD
      CALL PL(X0,Y0,XE,Y0,U)
      CALL PL(X0,Y0-3*SY,XE,Y0-3*SY,U)
      CALL PL(X0,Y0,X0,Y0-SY,U)
      CALL PL(X0,Y0-SY,XS,Y0-SY,U)
      CALL PL(XS,Y0,XS,Y0-SY,U)
      XM=XS+XR*SXD
      CALL PL(XM,Y0,XM,Y0+12,U)
      CALL CHAWTH(CWL)
      CALL HEIGHT(NINT(HL*U))
      CALL TXTJUS('CENT')
      CALL TXTJUS('MIDDLE')
      WRITE(LB,'(F3.1,A)') XR,'h'
      CALL MESSAG(LB,NINT(XM*U),NINT((Y0+DLB)*U))
      CALL DISFIN
      END
      SUBROUTINE PL(X1,Y1,X2,Y2,U)
      CALL LINE(NINT(X1*U),NINT(Y1*U),NINT(X2*U),NINT(Y2*U))
      RETURN
      END
