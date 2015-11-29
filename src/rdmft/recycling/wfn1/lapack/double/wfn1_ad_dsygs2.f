      SUBROUTINE WFN1_AD_DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      USE WFN1_AD1
      IMPLICIT NONE
#include "blas/double/intf_wfn1_ad_daxpy.fh"
#include "blas/double/intf_wfn1_ad_dscal.fh"
#include "blas/double/intf_wfn1_ad_dsyr2.fh"
#include "blas/double/intf_wfn1_ad_dtrmv.fh"
#include "blas/double/intf_wfn1_ad_dtrsv.fh"
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      TYPE(WFN1_AD_DBLE) :: A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYGS2 reduces a real symmetric-definite generalized eigenproblem
*  to standard form.
*
*  If ITYPE = 1, the problem is A*x = lambda*B*x,
*  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
*
*  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
*  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T *A*L.
*
*  B must have been previously factorized as U**T *U or L*L**T by DPOTRF.
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
*          = 2 or 3: compute U*A*U**T or L**T *A*L.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored, and how B has been factorized.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the transformed matrix, stored in the
*          same format as A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
*          The triangular factor from the Cholesky factorization of B,
*          as returned by DPOTRF.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      TYPE(WFN1_AD_DBLE) :: ONE
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K
      TYPE(WFN1_AD_DBLE) :: AKK, BKK, CT
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
c     INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      ONE  = 1.0d0
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGS2', -INFO )
         RETURN
      END IF
*
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
*
*           Compute inv(U**T)*A*inv(U)
*
            DO 10 K = 1, N
*
*              Update the upper triangle of A(k:n,k:n)
*
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2.0d0
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL WFN1_AD_DSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL WFN1_AD_DAXPY( N-K, CT, B( K, K+1 ), LDB,
     $                 A( K, K+1 ), LDA )
                  CALL WFN1_AD_DSYR2( UPLO, N-K, -ONE, A( K, K+1 ), LDA,
     $                 B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )
                  CALL WFN1_AD_DAXPY( N-K, CT, B( K, K+1 ), LDB,
     $                 A( K, K+1 ), LDA )
                  CALL WFN1_AD_DTRSV( UPLO, 'Transpose', 'Non-unit',
     $                 N-K, B( K+1, K+1 ), LDB, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
*
*           Compute inv(L)*A*inv(L**T)
*
            DO 20 K = 1, N
*
*              Update the lower triangle of A(k:n,k:n)
*
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2.0d0
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL WFN1_AD_DSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL WFN1_AD_DAXPY( N-K, CT, B( K+1, K ), 1,
     $                 A( K+1, K ), 1 )
                  CALL WFN1_AD_DSYR2( UPLO, N-K, -ONE, A( K+1, K ), 1,
     $                 B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL WFN1_AD_DAXPY( N-K, CT, B( K+1, K ), 1,
     $                 A( K+1, K ), 1 )
                  CALL WFN1_AD_DTRSV( UPLO, 'No transpose', 'Non-unit',
     $                 N-K, B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
*
*           Compute U*A*U**T
*
            DO 30 K = 1, N
*
*              Update the upper triangle of A(1:k,1:k)
*
               AKK = A( K, K )
               BKK = B( K, K )
               CALL WFN1_AD_DTRMV( UPLO, 'No transpose', 'Non-unit',
     $              K-1, B, LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL WFN1_AD_DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL WFN1_AD_DSYR2( UPLO, K-1, ONE, A( 1, K ), 1,
     $              B( 1, K ), 1, A, LDA )
               CALL WFN1_AD_DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL WFN1_AD_DSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2.0d0
   30       CONTINUE
         ELSE
*
*           Compute L**T *A*L
*
            DO 40 K = 1, N
*
*              Update the lower triangle of A(1:k,1:k)
*
               AKK = A( K, K )
               BKK = B( K, K )
               CALL WFN1_AD_DTRMV( UPLO, 'Transpose', 'Non-unit', K-1,
     $              B, LDB, A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL WFN1_AD_DAXPY( K-1, CT, B( K, 1 ), LDB,
     $              A( K, 1 ), LDA )
               CALL WFN1_AD_DSYR2( UPLO, K-1, ONE, A( K, 1 ), LDA,
     $              B( K, 1 ), LDB, A, LDA )
               CALL WFN1_AD_DAXPY( K-1, CT, B( K, 1 ), LDB,
     $              A( K, 1 ), LDA )
               CALL WFN1_AD_DSCAL( K-1, BKK, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2.0d0
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DSYGS2
*
      END
