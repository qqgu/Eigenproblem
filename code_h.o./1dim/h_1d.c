/********This is a code for one dimension harmonic oscillator eigen-problem************/
#include <stdio.h>
//lapacke headers
#include "mkl.h"

int main(){
    int matrix_layout = LAPACK_COL_MAJOR;
    char jobz = 'V';
    char uplo='U';
    int M,M2;
    double xmax,xmin;//(range of x)
    double h,R;
    M = 800;
    M2 = M/2;
    xmax = 10.0;
    xmin = -10.0;
    h = (xmax-xmin)/((double) M);
    R = 1./(h*h);
    int lda = M-1;
    int order = M-1;
    double *A,*w;
    FILE *s,*v;
    s=fopen("eigva.dat","w");
    v=fopen("eigve.dat","w");

    A = (double *)mkl_malloc( order*order*sizeof( double ), 64);
    w = (double *)mkl_malloc( order*sizeof( double ), 64);
    int i,j;
    for(i=0;i<order;i++)
    {
            A[i*order+i]=2.0*R+((i+1-M2)*h)*((i+1-M2)*h); //(1-M2) -> M2-1 *h
            if(i<order-1) A[i*order+i+1]=-R;
            if(i>0) A[i*order+i-1]=-R;
            w[i]=0;
    }

    int info = LAPACKE_dsyev(matrix_layout,jobz,uplo,order,A,lda,w);
    if(info==0){
        for(i=0;i<order;i++){
          //  fprintf(s,"eigenvalue %d:\n",i);
            fprintf(s,"%.6g\n",w[i]);       //print eigenvalue.
          // fprintf(v,"right eigenvector: \n");
            for(j=0;j<order;j++)
                fprintf(v,"%.6g \t",A[i*order+j]);  //print eigenvector.
            fprintf(v,"\n");
        }
        printf("SUCCESS\n");
        
    }
    fclose(s);
    fclose(v);
    mkl_free(A);
    mkl_free(w);

    return 0;
}
