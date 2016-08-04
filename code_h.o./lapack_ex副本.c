#include <stdio.h>

//lapacke headers
#include "mkl.h"

int main(){
    int matrix_order = LAPACK_COL_MAJOR;
    char jobvl = 'N';
    char jobvr = 'V';
    int n = 4;
    double A[16] = {
         0.35,  0.09, -0.44,  0.25,
         0.45,  0.07, -0.33, -0.32,
        -0.14, -0.54, -0.03, -0.13,
        -0.17,  0.35,  0.17,  0.11
        };
    int lda = n;
    double wr[4] = {0,0,0,0};
    double wi[4] = {0,0,0,0};
    double vl[16];
    int ldvl = 4;
    double vr[16];
    int ldvr = 4;

    int info = LAPACKE_dgeev(matrix_order,jobvl,jobvr,n,A,lda,wr,wi,vl,ldvl,vr,ldvr);
    if(info==0){
        int i = 0;
        int j = 0;
        for(i=0;i<n;i++){
            printf("eigenvalue %d:\n",i);
            printf("%.6g + i %.6g \t",wr[i],wi[i]);
            printf("right eigenvector: ");
            for(j=0;j<ldvr;j++)
                printf("%.6g \t",vr[i*4+j]);
            printf("\n");
        }
        printf("SUCCESS\n");
    }

    return 0;
}
