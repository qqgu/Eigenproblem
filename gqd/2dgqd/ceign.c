#include <stdio.h>
#include "mkl.h"

int main(){
    int matrix_layout = LAPACK_COL_MAJOR;
    char jobz = 'V';
    char uplo = 'U';
    int n = 4;
    int lda = n;
	MKL_Complex16  *A = (MKL_Complex16 *) mkl_malloc(sizeof(MKL_Complex16) * n * n,64);
    double *w;
    w=mkl_malloc(n * sizeof(double), 64);
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            A[i*n+j].real=i+j;
            if(i<j) A[i*n+j].imag= i+j;
            if(i>j) A[i*n+j].imag= -(i+j);
        }
    }

	 int info=LAPACKE_zheev(matrix_layout, jobz, uplo, n, A, lda, w);
    
	 if(info==0)
    {
        for(i=0;i<n;i++){
            printf("eigenvalue %d:\n",i);
            printf("%.6g \t",w[i]);
            printf("eigenvector: ");
            for(j=0;j<n;j++)
                printf("%.6g+%.6gi\t",A[i*n+j].real,A[i*n+j].imag);
            printf("\n");
        }
        printf("SUCCESS\n");
    }

    return 0;
}
