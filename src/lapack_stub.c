/*
 This code reproduces two of the functions of lapack with an implementation
 using gsl. This allows the use of these functions without the need to link to
 any fortran code. For compatibility, all of the same arguments are taken, but
 some of the arguments may not be used. For instance dposv_ takes an argument U
 or L but will actually return the full solution, not simply the upper or
 lower part.
*/
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

void dposv_(char *UPLO, long *N, long *NRHS, double *A, long *LDA, double *B, long *LDB, long *info)
{
        /*
          Note:
          Linear algebra function from lapack re-implemented with gsl, what
          follows is the argument description from the origional function

          DPOSV - compute the solution to a real system of linear equations A * X = B,

          ARGUMENTS
          UPLO    (input) CHARACTER*1
                  = 'U':  Upper triangle of A is stored;
                  = 'L':  Lower triangle of A is stored.

          N       (input) INTEGER
                  The number of linear equations, i.e., the order of
                  the matrix A.  N >= 0.

          NRHS    (input) INTEGER
                  The number of right hand sides, i.e., the number of
                  columns of the matrix B.  NRHS >= 0.

          A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
                  On entry, the symmetric matrix A.  If UPLO = 'U',
                  the leading N-by-N upper triangular part of A con-
                  tains the upper triangular part of the matrix A, and
                  the strictly lower triangular part of A is not
                  referenced.  If UPLO = 'L', the leading N-by-N lower
                  triangular part of A contains the lower triangular
                  part of the matrix A, and the strictly upper tri-
                  angular part of A is not referenced.

                  On exit, if INFO = 0, the factor U or L from the
                  Cholesky factorization A = U**T*U or A = L*L**T.

          LDA     (input) INTEGER
                  The leading dimension of the array A.  LDA >=
                  max(1,N).

          B       (input/output) DOUBLE PRECISION array, dimension
                  On entry, the N-by-NRHS right hand side matrix B.
                  On exit, if INFO = 0, the N-by-NRHS solution matrix
                  X.

          LDB     (input) INTEGER
                  The leading dimension of the array B.  LDB >=
                  max(1,N).

          INFO    (output) INTEGER
                  = 0:  successful exit
                  < 0:  if INFO = -i, the i-th argument had an illegal
                  value
                  > 0:  if INFO = i, the leading minor of order i of A
                  is not positive definite, so the factorization could
                  not be completed, and the solution has not been com-
                  puted.

         */
	gsl_matrix_view m = gsl_matrix_view_array(A, *LDA, *N);
	gsl_vector_view b = gsl_vector_view_array(B, *LDB);

	*info = gsl_linalg_cholesky_decomp(&m.matrix);
	if (*info != 0)
	{
		return;
	}

	*info = gsl_linalg_cholesky_svx(&m.matrix, &b.vector);
	if (*info != 0)
	{
		return;
	}

        int i, j;

        fprintf(stderr, "%ld, %ld\n", *LDA, *N);
        for (i = 0; i < *LDA; i++) {
            for (j = 0; j < *N; j++) {
                fprintf(stderr, "%lf ", A[i * *LDA + j]);
            }
            fprintf(stderr, "\n");
        }

	*info = 0;
	return;
}


void dtrtri_(char *UPLO, char *DIAG, long *N, double *A, long *LDA, long *info)
{
        /*
        Note:
        Linear algebra function from lapack re-implemented with gsl, what
        follows is the argument description from the origional function

        DTRTRI computes the inverse of a real upper or lower tri-
        angular matrix A.

        This is the Level 3 BLAS version of the algorithm.

        ARGUMENTS
        UPLO    (input) CHARACTER*1
                = 'U':  A is upper triangular;
                = 'L':  A is lower triangular.

        DIAG    (input) CHARACTER*1
                = 'N':  A is non-unit triangular;
                = 'U':  A is unit triangular.

        N       (input) INTEGER
                The order of the matrix A.  N >= 0.

        A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
                On entry, the triangular matrix A.  If UPLO = 'U',
                the leading N-by-N upper triangular part of the
                array A contains the upper triangular matrix, and
                the strictly lower triangular part of A is not
                referenced.  If UPLO = 'L', the leading N-by-N lower
                triangular part of the array A contains the lower
                triangular matrix, and the strictly upper triangular
                part of A is not referenced.  If DIAG = 'U', the
                diagonal elements of A are also not referenced and
                are assumed to be 1.  On exit, the (triangular)
                inverse of the original matrix, in the same storage
                format.

        LDA     (input) INTEGER
                The leading dimension of the array A.  LDA >=
                max(1,N).

        INFO    (output) INTEGER
                = 0: successful exit

                < 0: if INFO = -i, the i-th argument had an illegal
                value
                > 0: if INFO = i, A(i,i) is exactly zero.  The tri-
                angular matrix is singular and its inverse can not
                be computed.

          */
	CBLAS_TRANSPOSE_t transa;
	long start, end, i, j;
	gsl_matrix_view m = gsl_matrix_view_array(A, *LDA, *N);
	gsl_vector *temp = gsl_vector_calloc(*LDA**LDA);
	gsl_vector_view view;

	// Check if the upper or lower triangle is to be inverted: these blocks
	// define some conditions that will always be true in these circumstances,
	// such as the start or ending position of copying the inverted row into
	// origional matrix. We also set if the matrix should be inverted for the
	// solution or not.
	if (*UPLO == 'L')
	{
		start = end = 0;
		transa = CblasNoTrans;
	}
	else
	{
		start = 1;
		end = *LDA;
		transa = CblasTrans;
	}

	for (i = 0; i < *LDA; i++)
	{
		temp->data[i**LDA+i] = 1;
		view = gsl_vector_subvector(temp, i**LDA, *LDA);
		*info = gsl_blas_dtrsv(CblasLower, transa, CblasNonUnit, &m.matrix, &view.vector);
		if (*UPLO == 'L')
		{
			end = i;
		}
		else
		{
			start = i+1;
		}

		for (j = start; j < end; j++)
		{
			temp->data[i**LDA+j] = A[i**LDA+j];
		}
	}
	for (i = 0; i < *LDA**LDA; i++)
	{
		A[i] = temp->data[i];
	}

	gsl_vector_free(temp);
	return;
}
