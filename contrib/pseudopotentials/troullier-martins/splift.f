       SUBROUTINE SPLIFT (X,Y,YP,YPP,N,W,IERR,ISX,A1,B1,AN,BN)
C
       implicit double precision (a-h, o-z)

       PARAMETER (FOUR=4.D0)
CRAY      PARAMETER (FOUR=4.0)
C
C  NJTJ
C  ###  CRAY CONVERSIONS
C  ###    1)Comment out the implicit double precision.
C  ###    2)Switch double precision parameter
C  ###      to single precision parameter
C  ###  CRAY CONVERSIONS
C  NJTJ
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
  237 FORMAT(F5.1,F5.1,2I5)
C
C     ABSTRACT
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
C
C       --INPUT--
C
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
C         Y    - ARRAY OF ORDINATES OF DATA
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
C                SET TO ONE FOR FASTER EXECUTION.
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
C                OF THE SPLINE AT THE END POINTS (SEE YPP).
C                THE END CONDITION CONSTRAINTS ARE
C                        YPP(1) = A1*YPP(2) + B1
C                AND
C                        YPP(N) = AN*YPP(N-1) + BN
C                WHERE
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
C
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
C
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
C                FOLLOWING CONDITIONS.
C
C                A1=-0.5
C
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
C
C                AN=-0.5
C
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
C
C       --OUTPUT--
C
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
C         IERR - A STATUS CODE
C              --NORMAL CODE
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
C              --ABNORMAL CODES
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.
C
C       --WORK--
C
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
       DIMENSION X(N),Y(N),YP(N),YPP(N),W(N,3)
C
       IF (N.LT.4) THEN
         IERR = 2
         RETURN
       ENDIF
       NM1  = N-1
       NM2  = N-2
       IF (ISX.GT.0) GO TO 40
       DO 5 I=2,N
         IF (X(I)-X(I-1) .LE. 0) THEN
           IERR = 3
           RETURN
         ENDIF
 5     CONTINUE
C
C     DEFINE THE TRIDIAGONAL MATRIX
C
       W(1,3) = X(2)-X(1)
       DO 10 I=2,NM1
         W(I,2) = W(I-1,3)
         W(I,3) = X(I+1)-X(I)
 10      W(I,1) = 2*(W(I,2)+W(I,3))
       W(1,1) = FOUR
       W(1,3) =-4*A1
       W(N,1) = FOUR
       W(N,2) =-4*AN
C
C     L U DECOMPOSITION
C
       DO 30 I=2,N
         W(I-1,3) = W(I-1,3)/W(I-1,1)
 30    W(I,1) = W(I,1) - W(I,2)*W(I-1,3)
C
C     DEFINE *CONSTANT* VECTOR
C
 40   YPP(1) = 4*B1
      DOLD = (Y(2)-Y(1))/W(2,2)
      DO 50 I=2,NM2
        DNEW   = (Y(I+1) - Y(I))/W(I+1,2)
        YPP(I) = 6*(DNEW - DOLD)
        YP(I)  = DOLD
 50   DOLD = DNEW
      DNEW = (Y(N)-Y(N-1))/(X(N)-X(N-1))
      YPP(NM1) = 6*(DNEW - DOLD)
      YPP(N) = 4*BN
      YP(NM1)= DOLD
      YP(N) = DNEW
C
C     FORWARD SUBSTITUTION
C
      YPP(1) = YPP(1)/W(1,1)
      DO 60 I=2,N
 60   YPP(I) = (YPP(I) - W(I,2)*YPP(I-1))/W(I,1)
C
C     BACKWARD SUBSTITUTION
C
       DO 70 J=1,NM1
         I = N-J
   70 YPP(I) = YPP(I) - W(I,3)*YPP(I+1)
C
C     COMPUTE FIRST DERIVATIVES
C
      YP(1) = (Y(2)-Y(1))/(X(2)-X(1)) - (X(2)-X(1))*(2*YPP(1)
     1  + YPP(2))/6
      DO 80 I=2,NM1
 80   YP(I) = YP(I) + W(I,2)*(YPP(I-1) + 2*YPP(I))/6
      YP(N) = YP(N) + (X(N)-X(NM1))*(YPP(NM1) + 2*YPP(N))/6
C
      IERR = 1
      RETURN
      END
       SUBROUTINE SPLINT (X,Y,YPP,N,XI,YI,YPI,YPPI,NI,KERR)
       implicit double precision (a-h,o-z)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C
C         SPLINT EVALUATES A CUBIC SPLINE AND ITS FIRST AND SECOND
C         DERIVATIVES AT THE ABSCISSAS IN XI.  THE SPLINE (WHICH
C         IS DEFINED BY X, Y, AND YPP) MAY HAVE BEEN DETERMINED BY
C         SPLIFT OR SMOO OR ANY OTHER SPLINE FITTING ROUTINE THAT
C         PROVIDES SECOND DERIVATIVES.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YPP(N), XI(NI), YI(NI), YPI(NI), YPPI(NI)
C
C       --INPUT--
C
C         X   - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C         Y   - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C         YPP - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C         N   - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C               THE ARRAYS X, Y, AND YPP MUST BE DIMENSIONED AT LEAST N
C               N MUST BE GREATER THAN OR EQUAL TO 2.
C         XI  - THE ABSCISSA OR ARRAY OF ABSCISSAS (IN ARBITRARY ORDER)
C               AT WHICH THE SPLINE IS TO BE EVALUATED.
C               EACH XI(K) THAT LIES BETWEEN X(1) AND X(N) IS A CASE OF
C               INTERPOLATION.  EACH XI(K) THAT DOES NOT LIE BETWEEN
C               X(1) AND X(N) IS A CASE OF EXTRAPOLATION.  BOTH CASES
C               ARE ALLOWED.  SEE DESCRIPTION OF KERR.
C         NI  - THE NUMBER OF ABSCISSAS AT WHICH THE SPLINE IS TO BE
C               EVALUATED.  IF NI IS GREATER THAN 1, THEN XI, YI, YPI,
C               AND YPPI MUST BE ARRAYS DIMENSIONED AT LEAST NI.
C               NI MUST BE GREATER THAN OR EQUAL TO 1.
C
C       --OUTPUT--
C
C         YI  - ARRAY OF VALUES OF THE SPLINE (ORDINATES) AT XI.
C         YPI - ARRAY OF VALUES OF THE FIRST DERIVATIVE OF SPLINE AT XI
C         YPPI- ARRAY OF VALUES OF SECOND DERIVATIVES OF SPLINE AT XI.
C         KERR- A STATUS CODE
C             --NORMAL CODES
C                1 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI USING ONLY INTERPOLATION.
C                2 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI, BUT AT LEAST ONE EXTRAPOLATION WAS PERFORMED.
C             -- ABNORMAL CODE
C                3 MEANS THAT THE REQUESTED NUMBER OF EVALUATIONS, NI,
C                  WAS NOT POSITIVE.
C
       DIMENSION X(N),Y(N),YPP(N),XI(NI),YI(NI),YPI(NI),YPPI(NI)
C
C     CHECK INPUT
C
      IF (NI) 1,1,2
 1    CONTINUE
C    1 CALL ERRCHK(67,67HIN SPLINT,  THE REQUESTED NUMBER OF INTERPOLATI
C     1NS WAS NOT POSITIVE)
      KERR = 3
      RETURN
    2 KERR = 1
      NM1= N-1
C
C     K IS INDEX ON VALUE OF XI BEING WORKED ON.  XX IS THAT VALUE.
C     I IS CURRENT INDEX INTO X ARRAY.
C
       K  = 1
       XX = XI(1)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IL = 1
       IR = N
C
C     BISECTION SEARCH
C
   10 I  = (IL+IR)/2
       IF (I.EQ.IL) GO TO 100
       IF (XX-X(I)) 20,100,30
   20 IR = I
       GO TO 10
   30 IL = I
       GO TO 10
C
C     LINEAR FORWARD SEARCH
C
   50 IF (XX-X(I+1)) 100,100,60
   60 IF (I.GE.NM1) GO TO 80
       I  = I+1
       GO TO 50
C
C     EXTRAPOLATION
C
   80 KERR = 2
      I  = NM1
      GO TO 100
   90 KERR = 2
      I  = 1
C
C     INTERPOLATION
C
  100 H  = X(I+1) - X(I)
       H2 = H*H
       XR = (X(I+1)-XX)/H
       XR2= XR*XR
       XR3= XR*XR2
       XL = (XX-X(I))/H
       XL2= XL*XL
       XL3= XL*XL2
       YI(K) = Y(I)*XR + Y(I+1)*XL
     1       -H2*(YPP(I)*(XR-XR3) + YPP(I+1)*(XL-XL3))/6.0D0
       YPI(K) = (Y(I+1)-Y(I))/H
     1 +H*(YPP(I)*(1.0D0-3.0D0*XR2)-YPP(I+1)*(1.0D0-3.0D0*XL2))/6.0D0
       YPPI(K) = YPP(I)*XR + YPP(I+1)*XL
C
C     NEXT POINT
C
       IF (K.GE.NI) RETURN
       K = K+1
       XX = XI(K)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IF (XX-XI(K-1)) 110,100,50
  110 IL = 1
       IR = I+1
       GO TO 10
C
       END
