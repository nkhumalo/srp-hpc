      SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYGST reduces a real symmetric-definite generalized eigenproblem
*  to standard form.
*
*  If ITYPE = 1, the problem is A*x = lambda*B*x,
*  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
*
*  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
*  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
*
*  B must have been previously factorized as U**T*U or L*L**T by DPOTRF.
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
*          = 2 or 3: compute U*A*U**T or L**T*A*L.
*
*  UPLO    (input) CHARACTER
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored, and how B has been factorized.
*          = 'U':  Upper triangle of A is stored and B is factored as
*                  U**T*U;
*          = 'L':  Lower triangle of A is stored and B is factored as
*                  L*L**T.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
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
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSYGS2, DSYMM, DSYR2K, DTRMM, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
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
         CALL XERBLA( 'DSYGST', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DSYGST', UPLO, N, -1, -1, -1 )
*
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
*
*        Use blocked code
*
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
*
*              Compute inv(U')*A*inv(U)
*
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the upper triangle of A(k:n,k:n)
*
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Left', UPLO, 'Transpose', 'Non-unit',
     $                           KB, N-K-KB+1, ONE, B( K, K ), LDB,
     $                           A( K, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL DSYR2K( UPLO, 'Transpose', N-K-KB+1, KB, -ONE,
     $                            A( K, K+KB ), LDA, B( K, K+KB ), LDB,
     $                            ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL DTRSM( 'Right', UPLO, 'No transpose',
     $                           'Non-unit', KB, N-K-KB+1, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K, K+KB ),
     $                           LDA )
                  END IF
   10          CONTINUE
            ELSE
*
*              Compute inv(L)*A*inv(L')
*
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the lower triangle of A(k:n,k:n)
*
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                           N-K-KB+1, KB, ONE, B( K, K ), LDB,
     $                           A( K+KB, K ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL DSYR2K( UPLO, 'No transpose', N-K-KB+1, KB,
     $                            -ONE, A( K+KB, K ), LDA, B( K+KB, K ),
     $                            LDB, ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL DTRSM( 'Left', UPLO, 'No transpose',
     $                           'Non-unit', N-K-KB+1, KB, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K+KB, K ),
     $                           LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
*
*              Compute U*A*U'
*
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL DTRMM( 'Left', UPLO, 'No transpose', 'Non-unit',
     $                        K-1, KB, ONE, B, LDB, A( 1, K ), LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DSYR2K( UPLO, 'No transpose', K-1, KB, ONE,
     $                         A( 1, K ), LDA, B( 1, K ), LDB, ONE, A,
     $                         LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DTRMM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                        K-1, KB, ONE, B( K, K ), LDB, A( 1, K ),
     $                        LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
*
*              Compute L'*A*L
*
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL DTRMM( 'Right', UPLO, 'No transpose', 'Non-unit',
     $                        KB, K-1, ONE, B, LDB, A( K, 1 ), LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DSYR2K( UPLO, 'Transpose', K-1, KB, ONE,
     $                         A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A,
     $                         LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DTRMM( 'Left', UPLO, 'Transpose', 'Non-unit', KB,
     $                        K-1, ONE, B( K, K ), LDB, A( K, 1 ), LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
*
*     End of DSYGST
*
      END
