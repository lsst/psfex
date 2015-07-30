#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

//extern void dposv() asm ("dposv_");
//extern void dtrtri() asm ("dtrtri_");

void dposv_(char *UPLO, long *N, long *NRHS, double *A, long *LDA, double *B, long *LDB, long *info)
{
	//gsl_error_handler_t* newwerr = gsl_set_error_handler_off ();
	gsl_matrix_view m = gsl_matrix_view_array(A,*LDA,*N);
	gsl_vector_view b = gsl_vector_view_array(B,*LDB);

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

	*info = 0;
	return;
}


void dtrtri_(char *UPLO, char *DIAG, long * N, double * A, long * LDA, long * info)
{
	CBLAS_TRANSPOSE_t transa;
	long start,end,i,j;
	gsl_matrix_view m = gsl_matrix_view_array(A,*LDA,*N);
	gsl_vector * temp = gsl_vector_calloc(*LDA**LDA);
	gsl_vector_view view;

	//check if the upper or lower triangle is to be inverted
	//these blocks define some conditions that will always be true in
	//these circumstances. Such as the start or ending position of
	//copying the inverted row into origional matrix.
	//Finally we need to set if the matrix should be inverted for the
	//solution or not.
	if(*UPLO == 'L')
	{
		start  = 0;
		transa = CblasNoTrans;
	}
	else
	{
		end = *LDA;
		transa = CblasTrans;
	}

	for(i=0; i<*LDA;i++)
	{
		temp->data[i**LDA+i] = 1;
		view = gsl_vector_subvector(temp,i**LDA,*LDA);
		*info = gsl_blas_dtrsv(CblasLower,transa,CblasNonUnit,&m.matrix,
				       &view.vector);
		if(*UPLO == 'L')
		{
			end = i;
		}
		else
		{
			start = i+1;
		}
		for(j=0;j<*LDA;j++)
		{
			printf("%lf ",view.vector.data[j]);
		}

		for(j=start;j<end;j++)
		{
			printf("%ld ",j);
			temp->data[i**LDA+j] = A[i**LDA+j];
			printf("%lf ",temp->data[i**LDA+j]);
		}
		printf("\n");
	}
	for(i=0;i<*LDA**LDA;i++)
	{
		A[i] = temp->data[i];
	}

	gsl_vector_free(temp);
	return;
}


