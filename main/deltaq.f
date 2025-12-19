CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DELTAQ(IH,ICH,IPA,x,q,z,delta)
      IMPLICIT DOUBLE PRECISION (A-H,J-Z)
      double precision f(-6:6)
      COMMON /ORD/ IORD
      COMMON /POL/ IPOL
      COMMON /TAR/ ITAR 
      COMMON /PART/ IPART
C ***  CHOOSE TARGET : 0=proton; 1=neutron; 2=deuteron; 3=helium
      COMMON /HAD/ IHAD
      COMMON /CHAR/ ICHAR
      COMMON /KINVAR/ Q2,XB,ZHH
C **** IZ is set to 0 and ZMIN=1
      COMMON /ZINT/ IZ
      COMMON /CUT/ ZMIN
C******* COMMUNICATION WITH VEGAS **********************************
      COMMON /RESULT/ ERG1,ERG2,ERG3,ERG4
C*******************************************************************
      EXTERNAL SUBDELTA
      IORD=1
      IPOL=1
      ITAR=0
      IHAD=IH
      ICHAR=ICH
      IPART=IPA
      XB=x
      Q2=q
      ZHH=z
      IZ=0
      ZMIN = 1.D0
      CALL VEGAS(SUBDELTA,1.D-6,3,15000,3,0,0)
      delta = ERG1
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION FUNCTION SUBDELTA(XX)
C ********************************************************************
C ***  INTEGRAND IN THE (NEXT-TO-)LEADING ORDER CASE *****************
C ********************************************************************
       IMPLICIT DOUBLE PRECISION (A-H,J-Z)
       DIMENSION XX(3)
       COMMON /KINVAR/ Q2,XB,ZHH
       COMMON /PART/ IPART
       COMMON /HAD/ IHAD
       COMMON /CHAR/ ICHAR
       COMMON /ORD/ IORD
       COMMON /POL/ IPOL
       COMMON /TAR/ ITAR
       COMMON /ZINT/ IZ
       COMMON /CUT/ ZMIN

       PI = DACOS(-1.D0)
       CF = 4.D0/3.D0
       TR = 1.D0/2.D0
C********** JACOBIAN UND GRENZEN *************
C
C***** ZH - INTEGRATION *****
       IF (IZ.EQ.0) THEN
         ZH = ZHH
         J3 = 1.D0
       ELSE
         ZMAX = 0.9D0
         ZH = ZMIN+XX(3)*(ZMAX-ZMIN)
         J3 = ZMAX-ZMIN
       ENDIF
C
C***** X - INTEGRATION *****
       XDO = XB
       XUP = 1.D0
       X = XDO+XX(2)*(XUP-XDO)
       J2 = XUP-XDO
C
C***** Z - INTEGRATION *****
       ZDO = ZH
       ZUP = 1.D0
       Z = ZDO+XX(1)*(ZUP-ZDO)
       J1 = ZUP-ZDO
C
C ***  X - VALUES FOR PARTON DISTRIBUTIONS
       XP = XB/X
       XP0 = XB
       XLO = DLOG(1.D0-XDO)
       ZP = ZH/Z
       ZP0 = ZH
       ZLO = DLOG(1.D0-ZDO)
C
C ***  CALL OF PARTON DISTRIBUTIONS AND FRAGMENTATION FUNCTIONS :
       CALL PARCOMDEL (IORD,IPOL,ITAR,ICHAR,IPART,XP0,ZP0,Q2,QQ00,
     1     QG00,GQ00) 
       CALL PARCOMDEL (IORD,IPOL,ITAR,ICHAR,IPART,XP,ZP0,Q2,QQ0,
     1     QG0,GQ0)
       CALL PARCOMDEL (IORD,IPOL,ITAR,ICHAR,IPART,XP0,ZP,Q2,Q0Q,
     1      Q0G,G0Q)  
       CALL PARCOMDEL (IORD,IPOL,ITAR,ICHAR,IPART,XP,ZP,Q2,QQ,
     1     QG,GQ)
C     **********************************************************************
C ***  LEADING ORDER CONTRIB. :
         SUB1 = QQ00
C ***  DELTA(1-X)*DELTA(1-Z) CONTRIB. :
         SUB2 = CF * ( - 8.D0 + (XLO+ZLO)**2 ) * QQ00
C ***  DELTA(1-X) CONTRIB. :
         SUB3QQ = DLOG(1.D0-Z)/(1.D0-Z) * ( (1.D0+Z**2) * Q0Q/Z -
     1        2.D0 * QQ00 ) +
     2            (1.D0+Z**2)/(1.D0-Z)/Z * DLOG(Z) * Q0Q -
     3            (1.D0+Z)/Z * XLO * Q0Q + (1.D0-Z)/Z * Q0Q +
     4        2.D0/(1.D0-Z) * ( Q0Q/Z - QQ00 ) * XLO
         SUB3QG = ( (1.D0+(1.D0-Z)**2)/Z *
     1              ( DLOG( Z*(1.D0-Z) ) + XLO ) + Z )/Z * Q0G
         SUB3 = CF * ( SUB3QQ + SUB3QG )
C ***  DELTA(1-Z) CONTRIB. :
         SUB4QQ = DLOG(1.D0-X)/(1.D0-X) * ( (1.D0+X**2) * QQ0/X -
     1                                       2.D0 * QQ00 ) -
     2            (1.D0+X**2)/(1.D0-X)/X * DLOG(X) * QQ0 -
     3            (1.D0+X)/X * ZLO * QQ0 + (1.D0-X)/X * QQ0 +
     4             2.D0/(1.D0-X) * ( QQ0/X - QQ00 ) * ZLO
         IF (IPOL.EQ.0) THEN
            SUB4GQ = ( (X**2+(1.D0-X)**2) * ( DLOG((1.D0-X)/X)
     1                  + ZLO ) + 2.D0*X*(1.D0-X) )/X * GQ0
         ELSE
            SUB4GQ = ( (X**2-(1.D0-X)**2) * ( DLOG((1.D0-X)/X)
     1                  + ZLO ) + 2.D0*(1.D0-X) )/X * GQ0
         ENDIF
         SUB4 = CF * SUB4QQ + TR * SUB4GQ
C ***  REMAINING CONTRIBUTIONS :
         SUB5QQ = - (1.D0+Z)/Z/(1.D0-X) * ( QQ/X - Q0Q ) -
     1              (1.D0+X)/X/(1.D0-Z) * ( QQ/Z - QQ0 ) +
     2               2.D0 * (QQ/X/Z - Q0Q/Z - QQ0/X + QQ00)/
     3              (1.D0-X)/(1.D0-Z)
         IF (IPOL.EQ.0) THEN
            SUB5QQ = SUB5QQ + 2.D0 * (1.D0+X*Z) * QQ/X/Z
            SUB5QG = (1.D0+(1.D0-Z)**2)/Z/(1.D0-X) * (QG/X - Q0G)/Z +
     1               ( 2.D0 * (1.D0+X-X*Z) - (1.D0+X)/Z ) * QG/X/Z
            SUB5GQ = (X**2+(1.D0-X)**2) * ( (GQ/Z - GQ0)/(1.D0-Z)/X +
     1                                      (1.D0/Z-2.D0) * GQ/X/Z )
         ELSE
            SUB5QQ = SUB5QQ + 2.D0 * (X+Z) * QQ/X/Z
            SUB5QG = (1.D0+(1.D0-Z)**2)/Z/(1.D0-X) * (QG/X - Q0G)/Z +
     1               ( 2.D0 * (1.D0+X-Z) - (1.D0+X)/Z ) * QG/X/Z
            SUB5GQ = (X**2-(1.D0-X)**2) * ( (GQ/Z - GQ0)/(1.D0-Z)/X +
     1                                      (1.D0/Z-2.D0) * GQ/X/Z )
         ENDIF
         SUB5 = CF * (SUB5QQ + SUB5QG ) + TR * SUB5GQ
C
         ALPHAS=alphasPDF(SQRT(Q2))/(4*3.14159265359)

         FAC = DBLE(IORD) * ALPHAS/2.D0/PI
         SUBDELTA = SUB1 +
     1            FAC * ( SUB2 + J1 * SUB3 + J2 * SUB4 + J1*J2 *SUB5 )
         SUBDELTA = SUBDELTA * J3
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE PARCOMDEL (IORD,IPOL,ITAR,ICHAR,IPART,XP,ZP,Q2,QQ,
     1     QG,GQ) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character name*64
      double precision f(-6:6)
      character*20 lparm
      logical has_photon
      COMMON /HAD/ IHAD

      IS=2
      IO=IORD
      IC=ICHAR
      IH=IHAD
      IP=IPART

      call fDSS(IS,IH,IC,IO,ZP,Q2,DPU,DPUB,DPD,DPDB,DPS,DPSB,DPC,
     1     DPB,DPG)
      
      dpu=dpu/zp
      dpub=dpub/zp
      dpd=dpd/zp
      dpdb=dpdb/zp
      dps=dps/zp
      dpsb=dpsb/zp
      dpc=dpc/zp
      dpcb=dpc
      dpg=dpg/zp
C     
      IF (IPOL.EQ.0) THEN
         call evolvePDF(xp,DSQRT(Q2),f)
         
         UP =   f(2)
         DO =   f(1)
         UB =   f(-2)
         DB =   f(-1)
         ST =   f(3)
         STB =  f(-3)
         GL  =  f(0)
         CH  =  f(4)
         BO  =  f(5)
*     
      ELSE
         CALL DSSVFIT(XP,Q2,DUV,DDV,UB,DB,ST,GL,G1P,G1N)
C     MODE=1
C     CALL POLFIT(MODE,XP,Q2,DUV,DDV,UB,DB,ST,GL,
C     #              G1P,G1N)
         up=duv+ub
         do=ddv+db
      ENDIF
C     
C     *** PROTON :
      IF (IP.EQ.1) THEN
            QQP = (DO*DPD)/9.D0
            QGP = (DO*DPG)/9.D0
            GQP = 0
      ELSE IF (IP.EQ.-1) THEN
            QQP = (DB*DPDB)/9.D0
            QGP = (DB*DPG)/9.D0
            GQP = 0
      ELSE IF (IP.EQ.2) THEN
            QQP = (4.D0*UP*DPU)/9.D0
            QGP = (4.D0*UP*DPG)/9.D0
            GQP = 0
      ELSE IF (IP.EQ.-2) THEN
            QQP = (4.D0*UB*DPUB)/9.D0
            QGP = (4.D0*UB*DPG)/9.D0
            GQP = 0
      ELSE IF (IP.EQ.3) THEN
            QQP = (ST*DPS)/9.D0
            QGP = (ST*DPG)/9.D0
            GQP = 0
      ELSE IF (IP.EQ.-3) THEN
            QQP = (ST*DPSB)/9.D0
            QGP = (ST*DPG)/9.D0
            GQP = 0
      ELSE IF (IP.EQ.21) THEN
            QQP = 0
            QGP = 0
            GQP = GL*(4.D0 * (DPU+DPUB) + (DPD+DPDB) + (DPS+DPSB))/9.D0
      ENDIF   

C     *** NEUTRON :
C            QQN = ( 4.D0 * (DO*DPU+DB*DPUB) + (UP*DPD+UB*DPDB) +
C     1           (ST*DPS+ST*DPSB) )/9.D0
C            QGN = ( 4.D0 * (DO+DB) + (UP+UB) + (ST+ST) ) * DPG/9.D0
C            GQN = GQP
C     
      IF (ITAR.EQ.0) THEN
         QQ = QQP
         QG = QGP
         GQ = GQP
C            ELSE IF (ITAR.EQ.1) THEN
C               QQ = QQN
C               QG = QGN
C               GQ = GQN
C            ELSE IF (ITAR.EQ.2) THEN
C               IF (IPOL.EQ.0) THEN
C                  QQ = 0.5 * (QQP+QQN)
C                  QG = 0.5 * (QGP+QGN)
C                  GQ = 0.5 * (GQP+GQN)
C               ELSE
C                  QQ = 0.5 * (QQP+QQN) * (1.D0 - 1.5D0*0.058D0)
C                  QG = 0.5 * (QGP+QGN) * (1.D0 - 1.5D0*0.058D0)
C                  GQ = 0.5 * (GQP+GQN) * (1.D0 - 1.5D0*0.058D0)
C               ENDIF
C            ELSE
C               IF (IPOL.EQ.0) THEN
C                  QQ =  (2.0 * QQP+QQN) / 3.0
C                  QG =  (2.0 * QGP+QGN) / 3.0
C                  GQ =  (2.0 * GQP+GQN) / 3.0
C               ELSE
C                  QQ =  (2.0 * (-0.027D0) * QQP+ (0.865D0) * QQN) / 3.0
C                  QG =  (2.0 * (-0.027D0) * QGP+ (0.865D0) * QGN) / 3.0
C                  GQ =  (2.0 * (-0.027D0) * GQP+ (0.865D0) * GQN) / 3.0
C               ENDIF
      ENDIF
      QQ = QQ/XP
      QG = QG/XP
      GQ = GQ/XP
      RETURN
      END
