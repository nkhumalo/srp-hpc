*
* $Id$
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           sint.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Forward and backward, 1D real sin FFT
C
      SUBROUTINE SINT(N,X,WSAVE)
      DOUBLE PRECISION X(*),WSAVE(*)

      NP1 = N + 1
      IW1 = N/2 + 1
      IW2 = IW1 + NP1
      IW3 = IW2 + NP1
      CALL SINT1(N,X,WSAVE,WSAVE(IW1),WSAVE(IW2),WSAVE(IW3))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           sinti.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Initialization routine for SINT
C
      SUBROUTINE SINTI(N,WSAVE)
      DOUBLE PRECISION PI
      DOUBLE PRECISION PIMACH
      DOUBLE PRECISION DUM
      DOUBLE PRECISION DT
      DOUBLE PRECISION WSAVE(*)

      PI = PIMACH(DUM)
      IF (N.LE.1) RETURN
      NS2 = N/2
      NP1 = N + 1
      DT = PI/DBLE(NP1)
      DO 101 K = 1,NS2
          WSAVE(K) = 2.D0*SIN(K*DT)
  101 CONTINUE
      CALL RFFTI(NP1,WSAVE(NS2+1))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           sint1.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE SINT1(N,WAR,WAS,XH,X,IFAC)
      DOUBLE PRECISION SQRT3
      DOUBLE PRECISION XHOLD
      DOUBLE PRECISION T1
      DOUBLE PRECISION T2
      DOUBLE PRECISION WAR(*),WAS(*),X(*),XH(*)
      integer*4 n4
      INTEGER IFAC(*)
      DATA SQRT3/1.73205080756888D0/
C
C FFTPACK 5.0 auxiliary routine
C
      DO 100 I = 1,N
          XH(I) = WAR(I)
          WAR(I) = X(I)
  100 CONTINUE
#if 1
      n4=n-2
      IF (n4) 101,102,103
#else
      IF (N-2) 101,102,103
#endif
  101 XH(1) = XH(1) + XH(1)
      GO TO 106
  102 XHOLD = SQRT3* (XH(1)+XH(2))
      XH(2) = SQRT3* (XH(1)-XH(2))
      XH(1) = XHOLD
      GO TO 106
  103 NP1 = N + 1
      NS2 = N/2
      X(1) = 0.D0
      DO 104 K = 1,NS2
          KC = NP1 - K
          T1 = XH(K) - XH(KC)
          T2 = WAS(K)* (XH(K)+XH(KC))
          X(K+1) = T1 + T2
          X(KC+1) = T2 - T1
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN.NE.0) X(NS2+2) = 4.D0*XH(NS2+1)
      CALL RFFTF1(NP1,X,XH,WAR,IFAC)
      XH(1) = .5D0*X(1)
      DO 105 I = 3,N,2
          XH(I-1) = -X(I)
          XH(I) = XH(I-2) + X(I-1)
  105 CONTINUE
      IF (MODN.NE.0) GO TO 106
      XH(N) = -X(N+1)
  106 DO 107 I = 1,N
          X(I) = WAR(I)
          WAR(I) = XH(I)
  107 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           rffti.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Initialization routine for RFFTB, RFFTF
C
      SUBROUTINE RFFTI(N,WSAVE)
      DOUBLE PRECISION WSAVE(*)

      IF (N.EQ.1) RETURN
      CALL RFFTI1(N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           rfftf1.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE RFFTF1(N,C,CH,WA,IFAC)
      DOUBLE PRECISION CH(*),C(*),WA(*)
      INTEGER IFAC(*)
C
C FFTPACK 5.0 auxiliary routine
C
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1 = 1,NF
          KH = NF - K1
          IP = IFAC(KH+3)
          L1 = L2/IP
          IDO = N/L2
          IDL1 = IDO*L1
          IW = IW - (IP-1)*IDO
          NA = 1 - NA
          IF (IP.NE.4) GO TO 102
          IX2 = IW + IDO
          IX3 = IX2 + IDO
          IF (NA.NE.0) GO TO 101
          CALL RADF4(IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
          GO TO 110
  101     CALL RADF4(IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
          GO TO 110
  102     IF (IP.NE.2) GO TO 104
          IF (NA.NE.0) GO TO 103
          CALL RADF2(IDO,L1,C,CH,WA(IW))
          GO TO 110
  103     CALL RADF2(IDO,L1,CH,C,WA(IW))
          GO TO 110
  104     IF (IP.NE.3) GO TO 106
          IX2 = IW + IDO
          IF (NA.NE.0) GO TO 105
          CALL RADF3(IDO,L1,C,CH,WA(IW),WA(IX2))
          GO TO 110
  105     CALL RADF3(IDO,L1,CH,C,WA(IW),WA(IX2))
          GO TO 110
  106     IF (IP.NE.5) GO TO 108
          IX2 = IW + IDO
          IX3 = IX2 + IDO
          IX4 = IX3 + IDO
          IF (NA.NE.0) GO TO 107
          CALL RADF5(IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          GO TO 110
  107     CALL RADF5(IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          GO TO 110
  108     IF (IDO.EQ.1) NA = 1 - NA
          IF (NA.NE.0) GO TO 109
          CALL RADFG(IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
          NA = 1
          GO TO 110
  109     CALL RADFG(IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
          NA = 0
  110     L2 = L1
  111 CONTINUE
      IF (NA.EQ.1) RETURN
      DO 112 I = 1,N
          C(I) = CH(I)
  112 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           rffti1.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE RFFTI1(N,WA,IFAC)
      DOUBLE PRECISION TPI
      DOUBLE PRECISION PIMACH
      DOUBLE PRECISION DUM
      DOUBLE PRECISION ARGH
      DOUBLE PRECISION ARGLD
      DOUBLE PRECISION FI
      DOUBLE PRECISION ARG
      DOUBLE PRECISION WA(*)
      INTEGER IFAC(*),NTRYH(4)
      integer*4 n4
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C
C FFTPACK 5.0 auxiliary routine
C
      NL = N
      NF = 0
      NTRY = 0
      J = 0
  101 J = J + 1
#if 1
      n4=j-4
      IF (n4) 102,102,103
#else
      IF (J-4) 102,102,103
#endif
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY + 2
  104 NQ = NL/NTRY
      NR = NL - NTRY*NQ
#if 1
      n4=NR
      IF (n4) 101,105,101
#else
      IF (NR) 101,105,101
#endif
  105 NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY.NE.2) GO TO 107
      IF (NF.EQ.1) GO TO 107
      DO 106 I = 2,NF
          IB = NF - I + 2
          IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL.NE.1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 2.0D0*PIMACH(DUM)
      ARGH = TPI/DBLE(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1.EQ.0) RETURN
      DO 110 K1 = 1,NFM1
          IP = IFAC(K1+2)
          LD = 0
          L2 = L1*IP
          IDO = N/L2
          IPM = IP - 1
          DO 109 J = 1,IPM
              LD = LD + L1
              I = IS
              ARGLD = DBLE(LD)*ARGH
              FI = 0.D0
              DO 108 II = 3,IDO,2
                  I = I + 2
                  FI = FI + 1.D0
                  ARG = FI*ARGLD
                  WA(I-1) = COS(ARG)
                  WA(I) = SIN(ARG)
  108         CONTINUE
              IS = IS + IDO
  109     CONTINUE
          L1 = L2
  110 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           radf2.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE RADF2(IDO,L1,CC,CH,WA1)
      DOUBLE PRECISION TR2
      DOUBLE PRECISION TI2
      DOUBLE PRECISION CH(IDO,2,L1),CC(IDO,L1,2),WA1(*)
      integer*4 n4
C
C FFTPACK 5.0 auxiliary routine
C
      DO 101 K = 1,L1
          CH(1,1,K) = CC(1,K,1) + CC(1,K,2)
          CH(IDO,2,K) = CC(1,K,1) - CC(1,K,2)
  101 CONTINUE
#if 1
      n4=ido-2
      IF (n4) 107,105,102
#else
      IF (IDO-2) 107,105,102
#endif
  102 IDP2 = IDO + 2
      DO 104 K = 1,L1
          DO 103 I = 3,IDO,2
              IC = IDP2 - I
              TR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
              TI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
              CH(I,1,K) = CC(I,K,1) + TI2
              CH(IC,2,K) = TI2 - CC(I,K,1)
              CH(I-1,1,K) = CC(I-1,K,1) + TR2
              CH(IC-1,2,K) = CC(I-1,K,1) - TR2
  103     CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2).EQ.1) RETURN
  105 DO 106 K = 1,L1
          CH(1,2,K) = -CC(IDO,K,2)
          CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           radf3.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE RADF3(IDO,L1,CC,CH,WA1,WA2)
      DOUBLE PRECISION TAUR
      DOUBLE PRECISION TAUI
      DOUBLE PRECISION CR2
      DOUBLE PRECISION DR2
      DOUBLE PRECISION DI2
      DOUBLE PRECISION DR3
      DOUBLE PRECISION DI3
      DOUBLE PRECISION CI2
      DOUBLE PRECISION TR2
      DOUBLE PRECISION TI2
      DOUBLE PRECISION TR3
      DOUBLE PRECISION TI3
      DOUBLE PRECISION CH(IDO,3,L1),CC(IDO,L1,3),WA1(*),WA2(*)
      DATA TAUR,TAUI/-.5D0,.866025403784439D0/
C
C FFTPACK 5.0 auxiliary routine
C
      DO 101 K = 1,L1
          CR2 = CC(1,K,2) + CC(1,K,3)
          CH(1,1,K) = CC(1,K,1) + CR2
          CH(1,3,K) = TAUI* (CC(1,K,3)-CC(1,K,2))
          CH(IDO,2,K) = CC(1,K,1) + TAUR*CR2
  101 CONTINUE
      IF (IDO.EQ.1) RETURN
      IDP2 = IDO + 2
      DO 103 K = 1,L1
          DO 102 I = 3,IDO,2
              IC = IDP2 - I
              DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
              DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
              DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
              DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
              CR2 = DR2 + DR3
              CI2 = DI2 + DI3
              CH(I-1,1,K) = CC(I-1,K,1) + CR2
              CH(I,1,K) = CC(I,K,1) + CI2
              TR2 = CC(I-1,K,1) + TAUR*CR2
              TI2 = CC(I,K,1) + TAUR*CI2
              TR3 = TAUI* (DI2-DI3)
              TI3 = TAUI* (DR3-DR2)
              CH(I-1,3,K) = TR2 + TR3
              CH(IC-1,2,K) = TR2 - TR3
              CH(I,3,K) = TI2 + TI3
              CH(IC,2,K) = TI3 - TI2
  102     CONTINUE
  103 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           radf4.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE RADF4(IDO,L1,CC,CH,WA1,WA2,WA3)
      DOUBLE PRECISION HSQT2
      DOUBLE PRECISION TR1
      DOUBLE PRECISION TR2
      DOUBLE PRECISION CR2
      DOUBLE PRECISION CI2
      DOUBLE PRECISION CR3
      DOUBLE PRECISION CI3
      DOUBLE PRECISION CR4
      DOUBLE PRECISION CI4
      DOUBLE PRECISION TR4
      DOUBLE PRECISION TI1
      DOUBLE PRECISION TI4
      DOUBLE PRECISION TI2
      DOUBLE PRECISION TI3
      DOUBLE PRECISION TR3
      DOUBLE PRECISION CC(IDO,L1,4),CH(IDO,4,L1),WA1(*),WA2(*),WA3(*)
      integer*4 n4
      DATA HSQT2/.7071067811865475D0/
C
C FFTPACK 5.0 auxiliary routine
C
      DO 101 K = 1,L1
          TR1 = CC(1,K,2) + CC(1,K,4)
          TR2 = CC(1,K,1) + CC(1,K,3)
          CH(1,1,K) = TR1 + TR2
          CH(IDO,4,K) = TR2 - TR1
          CH(IDO,2,K) = CC(1,K,1) - CC(1,K,3)
          CH(1,3,K) = CC(1,K,4) - CC(1,K,2)
  101 CONTINUE
#if 1
      n4=ido-2
      IF (n4) 107,105,102
#else
      IF (IDO-2) 107,105,102
#endif
  102 IDP2 = IDO + 2
      DO 104 K = 1,L1
          DO 103 I = 3,IDO,2
              IC = IDP2 - I
              CR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
              CI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
              CR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
              CI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
              CR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
              CI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
              TR1 = CR2 + CR4
              TR4 = CR4 - CR2
              TI1 = CI2 + CI4
              TI4 = CI2 - CI4
              TI2 = CC(I,K,1) + CI3
              TI3 = CC(I,K,1) - CI3
              TR2 = CC(I-1,K,1) + CR3
              TR3 = CC(I-1,K,1) - CR3
              CH(I-1,1,K) = TR1 + TR2
              CH(IC-1,4,K) = TR2 - TR1
              CH(I,1,K) = TI1 + TI2
              CH(IC,4,K) = TI1 - TI2
              CH(I-1,3,K) = TI4 + TR3
              CH(IC-1,2,K) = TR3 - TI4
              CH(I,3,K) = TR4 + TI3
              CH(IC,2,K) = TR4 - TI3
  103     CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2).EQ.1) RETURN
  105 CONTINUE
      DO 106 K = 1,L1
          TI1 = -HSQT2* (CC(IDO,K,2)+CC(IDO,K,4))
          TR1 = HSQT2* (CC(IDO,K,2)-CC(IDO,K,4))
          CH(IDO,1,K) = TR1 + CC(IDO,K,1)
          CH(IDO,3,K) = CC(IDO,K,1) - TR1
          CH(1,2,K) = TI1 - CC(IDO,K,3)
          CH(1,4,K) = TI1 + CC(IDO,K,3)
  106 CONTINUE
  107 RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           radf5.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE RADF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      DOUBLE PRECISION TR11
      DOUBLE PRECISION TI11
      DOUBLE PRECISION TR12
      DOUBLE PRECISION TI12
      DOUBLE PRECISION CR2
      DOUBLE PRECISION CI5
      DOUBLE PRECISION CR3
      DOUBLE PRECISION CI4
      DOUBLE PRECISION DR2
      DOUBLE PRECISION DI2
      DOUBLE PRECISION DR3
      DOUBLE PRECISION DI3
      DOUBLE PRECISION DR4
      DOUBLE PRECISION DI4
      DOUBLE PRECISION DR5
      DOUBLE PRECISION DI5
      DOUBLE PRECISION CR5
      DOUBLE PRECISION CI2
      DOUBLE PRECISION CR4
      DOUBLE PRECISION CI3
      DOUBLE PRECISION TR2
      DOUBLE PRECISION TI2
      DOUBLE PRECISION TR3
      DOUBLE PRECISION TI3
      DOUBLE PRECISION TR5
      DOUBLE PRECISION TI5
      DOUBLE PRECISION TR4
      DOUBLE PRECISION TI4
      DOUBLE PRECISION CC(IDO,L1,5),CH(IDO,5,L1),WA1(*),WA2(*),WA3(*),
     +                 WA4(*)
      DATA TR11,TI11,TR12,TI12/.309016994374947D0,.951056516295154D0,
     +     -.809016994374947D0,.587785252292473D0/
C
C FFTPACK 5.0 auxiliary routine
C
      DO 101 K = 1,L1
          CR2 = CC(1,K,5) + CC(1,K,2)
          CI5 = CC(1,K,5) - CC(1,K,2)
          CR3 = CC(1,K,4) + CC(1,K,3)
          CI4 = CC(1,K,4) - CC(1,K,3)
          CH(1,1,K) = CC(1,K,1) + CR2 + CR3
          CH(IDO,2,K) = CC(1,K,1) + TR11*CR2 + TR12*CR3
          CH(1,3,K) = TI11*CI5 + TI12*CI4
          CH(IDO,4,K) = CC(1,K,1) + TR12*CR2 + TR11*CR3
          CH(1,5,K) = TI12*CI5 - TI11*CI4
  101 CONTINUE
      IF (IDO.EQ.1) RETURN
      IDP2 = IDO + 2
      DO 103 K = 1,L1
          DO 102 I = 3,IDO,2
              IC = IDP2 - I
              DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
              DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
              DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
              DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
              DR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
              DI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
              DR5 = WA4(I-2)*CC(I-1,K,5) + WA4(I-1)*CC(I,K,5)
              DI5 = WA4(I-2)*CC(I,K,5) - WA4(I-1)*CC(I-1,K,5)
              CR2 = DR2 + DR5
              CI5 = DR5 - DR2
              CR5 = DI2 - DI5
              CI2 = DI2 + DI5
              CR3 = DR3 + DR4
              CI4 = DR4 - DR3
              CR4 = DI3 - DI4
              CI3 = DI3 + DI4
              CH(I-1,1,K) = CC(I-1,K,1) + CR2 + CR3
              CH(I,1,K) = CC(I,K,1) + CI2 + CI3
              TR2 = CC(I-1,K,1) + TR11*CR2 + TR12*CR3
              TI2 = CC(I,K,1) + TR11*CI2 + TR12*CI3
              TR3 = CC(I-1,K,1) + TR12*CR2 + TR11*CR3
              TI3 = CC(I,K,1) + TR12*CI2 + TR11*CI3
              TR5 = TI11*CR5 + TI12*CR4
              TI5 = TI11*CI5 + TI12*CI4
              TR4 = TI12*CR5 - TI11*CR4
              TI4 = TI12*CI5 - TI11*CI4
              CH(I-1,3,K) = TR2 + TR5
              CH(IC-1,2,K) = TR2 - TR5
              CH(I,3,K) = TI2 + TI5
              CH(IC,2,K) = TI5 - TI2
              CH(I-1,5,K) = TR3 + TR4
              CH(IC-1,4,K) = TR3 - TR4
              CH(I,5,K) = TI3 + TI4
              CH(IC,4,K) = TI4 - TI3
  102     CONTINUE
  103 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      File:           radfg.f
C
C      Library:        FFTPACK 4.1
C
C      Author:         Paul N. Swarztrauber
C                      National Center for Atmospheric Research
C                      PO 3000, Boulder, Colorado
C
C      Date:           Wed Mar 29 18:31:13 MST 1995
C
C      Description:    Lower-level auxiliary routine
C
      SUBROUTINE RADFG(IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      DOUBLE PRECISION TPI
      DOUBLE PRECISION PIMACH
      DOUBLE PRECISION DUM
      DOUBLE PRECISION ARG
      DOUBLE PRECISION DCP
      DOUBLE PRECISION DSP
      DOUBLE PRECISION AR1
      DOUBLE PRECISION AI1
      DOUBLE PRECISION AR1H
      DOUBLE PRECISION DC2
      DOUBLE PRECISION DS2
      DOUBLE PRECISION AR2
      DOUBLE PRECISION AI2
      DOUBLE PRECISION AR2H
      DOUBLE PRECISION CH(IDO,L1,IP),CC(IDO,IP,L1),C1(IDO,L1,IP),
     +                 C2(IDL1,IP),CH2(IDL1,IP),WA(*)
C
C FFTPACK 5.0 auxiliary routine
C
      TPI = 2.0D0*PIMACH(DUM)
      ARG = TPI/DBLE(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP + 2
      IDP2 = IDO + 2
      NBD = (IDO-1)/2
      IF (IDO.EQ.1) GO TO 119
      DO 101 IK = 1,IDL1
          CH2(IK,1) = C2(IK,1)
  101 CONTINUE
      DO 103 J = 2,IP
          DO 102 K = 1,L1
              CH(1,K,J) = C1(1,K,J)
  102     CONTINUE
  103 CONTINUE
      IF (NBD.GT.L1) GO TO 107
      IS = -IDO
      DO 106 J = 2,IP
          IS = IS + IDO
          IDIJ = IS
          DO 105 I = 3,IDO,2
              IDIJ = IDIJ + 2
              DO 104 K = 1,L1
                  CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J) +
     +                          WA(IDIJ)*C1(I,K,J)
                  CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J) -
     +                        WA(IDIJ)*C1(I-1,K,J)
  104         CONTINUE
  105     CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J = 2,IP
          IS = IS + IDO
          DO 109 K = 1,L1
              IDIJ = IS
              DO 108 I = 3,IDO,2
                  IDIJ = IDIJ + 2
                  CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J) +
     +                          WA(IDIJ)*C1(I,K,J)
                  CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J) -
     +                        WA(IDIJ)*C1(I-1,K,J)
  108         CONTINUE
  109     CONTINUE
  110 CONTINUE
  111 IF (NBD.LT.L1) GO TO 115
      DO 114 J = 2,IPPH
          JC = IPP2 - J
          DO 113 K = 1,L1
              DO 112 I = 3,IDO,2
                  C1(I-1,K,J) = CH(I-1,K,J) + CH(I-1,K,JC)
                  C1(I-1,K,JC) = CH(I,K,J) - CH(I,K,JC)
                  C1(I,K,J) = CH(I,K,J) + CH(I,K,JC)
                  C1(I,K,JC) = CH(I-1,K,JC) - CH(I-1,K,J)
  112         CONTINUE
  113     CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J = 2,IPPH
          JC = IPP2 - J
          DO 117 I = 3,IDO,2
              DO 116 K = 1,L1
                  C1(I-1,K,J) = CH(I-1,K,J) + CH(I-1,K,JC)
                  C1(I-1,K,JC) = CH(I,K,J) - CH(I,K,JC)
                  C1(I,K,J) = CH(I,K,J) + CH(I,K,JC)
                  C1(I,K,JC) = CH(I-1,K,JC) - CH(I-1,K,J)
  116         CONTINUE
  117     CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK = 1,IDL1
          C2(IK,1) = CH2(IK,1)
  120 CONTINUE
  121 DO 123 J = 2,IPPH
          JC = IPP2 - J
          DO 122 K = 1,L1
              C1(1,K,J) = CH(1,K,J) + CH(1,K,JC)
              C1(1,K,JC) = CH(1,K,JC) - CH(1,K,J)
  122     CONTINUE
  123 CONTINUE
C
      AR1 = 1.D0
      AI1 = 0.D0
      DO 127 L = 2,IPPH
          LC = IPP2 - L
          AR1H = DCP*AR1 - DSP*AI1
          AI1 = DCP*AI1 + DSP*AR1
          AR1 = AR1H
          DO 124 IK = 1,IDL1
              CH2(IK,L) = C2(IK,1) + AR1*C2(IK,2)
              CH2(IK,LC) = AI1*C2(IK,IP)
  124     CONTINUE
          DC2 = AR1
          DS2 = AI1
          AR2 = AR1
          AI2 = AI1
          DO 126 J = 3,IPPH
              JC = IPP2 - J
              AR2H = DC2*AR2 - DS2*AI2
              AI2 = DC2*AI2 + DS2*AR2
              AR2 = AR2H
              DO 125 IK = 1,IDL1
                  CH2(IK,L) = CH2(IK,L) + AR2*C2(IK,J)
                  CH2(IK,LC) = CH2(IK,LC) + AI2*C2(IK,JC)
  125         CONTINUE
  126     CONTINUE
  127 CONTINUE
      DO 129 J = 2,IPPH
          DO 128 IK = 1,IDL1
              CH2(IK,1) = CH2(IK,1) + C2(IK,J)
  128     CONTINUE
  129 CONTINUE
C
      IF (IDO.LT.L1) GO TO 132
      DO 131 K = 1,L1
          DO 130 I = 1,IDO
              CC(I,1,K) = CH(I,K,1)
  130     CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I = 1,IDO
          DO 133 K = 1,L1
              CC(I,1,K) = CH(I,K,1)
  133     CONTINUE
  134 CONTINUE
  135 DO 137 J = 2,IPPH
          JC = IPP2 - J
          J2 = J + J
          DO 136 K = 1,L1
              CC(IDO,J2-2,K) = CH(1,K,J)
              CC(1,J2-1,K) = CH(1,K,JC)
  136     CONTINUE
  137 CONTINUE
      IF (IDO.EQ.1) RETURN
      IF (NBD.LT.L1) GO TO 141
      DO 140 J = 2,IPPH
          JC = IPP2 - J
          J2 = J + J
          DO 139 K = 1,L1
              DO 138 I = 3,IDO,2
                  IC = IDP2 - I
                  CC(I-1,J2-1,K) = CH(I-1,K,J) + CH(I-1,K,JC)
                  CC(IC-1,J2-2,K) = CH(I-1,K,J) - CH(I-1,K,JC)
                  CC(I,J2-1,K) = CH(I,K,J) + CH(I,K,JC)
                  CC(IC,J2-2,K) = CH(I,K,JC) - CH(I,K,J)
  138         CONTINUE
  139     CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J = 2,IPPH
          JC = IPP2 - J
          J2 = J + J
          DO 143 I = 3,IDO,2
              IC = IDP2 - I
              DO 142 K = 1,L1
                  CC(I-1,J2-1,K) = CH(I-1,K,J) + CH(I-1,K,JC)
                  CC(IC-1,J2-2,K) = CH(I-1,K,J) - CH(I-1,K,JC)
                  CC(I,J2-1,K) = CH(I,K,J) + CH(I,K,JC)
                  CC(IC,J2-2,K) = CH(I,K,JC) - CH(I,K,J)
  142         CONTINUE
  143     CONTINUE
  144 CONTINUE
      RETURN
      END
