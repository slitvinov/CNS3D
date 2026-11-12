      SUBROUTINE PATT(IP,NP,U)
      INTEGER IP(*),JP(8)
      IF(NP.EQ.0) THEN
        CALL SOLID
        RETURN
      ENDIF
      DO 10 I=1,NP
  10  JP(I)=NINT(IP(I)*U)
      CALL MYLINE(JP,NP)
      RETURN
      END
      SUBROUTINE MARK2(NX,NY,R1,R2,R3,U)
      CALL LINWID(1)
      CALL CIRCSP(1)
      CALL SHDPAT(16)
      CALL CIRCLE(NX,NY,NINT(R1*U))
      CALL COLOR('BACK')
      CALL CIRCLE(NX,NY,NINT(R2*U))
      CALL COLOR('FORE')
      CALL CIRCLE(NX,NY,NINT(R3*U))
      CALL SHDPAT(0)
      RETURN
      END
      SUBROUTINE MYLEG(LAB,NLAB,KIND,IP,NP,XS,SLEN,TOFF,Y1,GLIN,GENT,
     1                 SDY,TDY,R1,R2,R3,SXO,SYO,WLIN,U)
      CHARACTER*(*) LAB(NLAB)
      INTEGER KIND(NLAB),IP(8,NLAB),NP(NLAB)
      NL=0
      DO 20 K=1,NLAB
      Y=Y1+NL*GLIN+(K-1)*GENT
      IF(KIND(K).EQ.1) THEN
        CALL LINWID(NINT(WLIN*U))
        CALL PATT(IP(1,K),NP(K),U)
        CALL LINE(NINT(XS*U),NINT((Y+SDY)*U),NINT((XS+SLEN)*U),
     1            NINT((Y+SDY)*U))
        CALL SOLID
      ELSE IF(KIND(K).EQ.2) THEN
        CALL MARK2(NINT((XS+SLEN/2+SXO)*U),NINT((Y+SDY+SYO)*U),
     1             R1,R2,R3,U)
      ENDIF
      I0=1
      DO 10 I=1,LEN_TRIM(LAB(K))+1
      IF(I.GT.LEN_TRIM(LAB(K)).OR.LAB(K)(I:I).EQ.'|') THEN
        CALL MESSAG(LAB(K)(I0:I-1),NINT((XS+SLEN+TOFF)*U),
     1              NINT((Y+TDY)*U))
        IF(I.LE.LEN_TRIM(LAB(K))) THEN
          Y=Y+GLIN
          NL=NL+1
        ENDIF
        I0=I+1
      ENDIF
  10  CONTINUE
  20  CONTINUE
      RETURN
      END
      SUBROUTINE RDXY(FILE,X,Y,N)
      CHARACTER*(*) FILE
      REAL X(*),Y(*)
      OPEN(1,FILE=FILE)
      READ(1,*) N
      DO 10 I=1,N
  10  READ(1,*) X(I),Y(I)
      CLOSE(1)
      RETURN
      END
