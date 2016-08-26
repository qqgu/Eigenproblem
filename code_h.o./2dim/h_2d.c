/********This is a code for one dimension harmonic oscillator eigen-problem************/
#include <stdio.h>
//lapacke headers
#include "mkl.h"

int main(){
    int matrix_layout = LAPACK_ROW_MAJOR;
    char jobz = 'V';
    char uplo='U';
    int M,M2,N,N2;
    double xmax,xmin,ymax,ymin;//(range of space)
    double hx,hy,Rx,Ry;
    M = 50;
    N = 50;
    M2 = M/2;
    N2 = N/2;
    xmax = 10.0;
    xmin = -10.0;
    ymax = 10.0;
    ymin = -10.0;
    
    hx = (xmax-xmin)/((double) N);
    hy = (ymax-ymin)/((double) M);
    Rx = 1.0/(hx*hx);
    Ry = 1.0/(hy*hy);
    int lda = (M-1)*(N-1);
    int order = (M-1)*(N-1);
    double *A,*w,*D,*C;
    FILE *s,*v;
    s=fopen("eigva.dat","w");
    v=fopen("eigve.dat","w");

    A = (double *)mkl_malloc( order*order*sizeof( double ), 64);
    w = (double *)mkl_malloc( order*sizeof( double ), 64);
    D = (double *)mkl_malloc( (M-1)*(M-1)*sizeof( double ), 64);
    C = (double *)mkl_malloc( (M-1)*(M-1)*sizeof( double ), 64);
    
    int i,j,m,n,mx,my,nx,ny;
    int di,dj;
    for (m=0;m<M-1;m++){
        D[m+m*(M-1)]=-Ry;
    }
    for(m=0;m<N-2;m++){
        mx = (M-1)*m;
        my = (M-1)*m+M-1;
        nx = (M-1)*m+M-1;
        ny = (M-1)*m+M-1+M-1;
        for (i=mx;i<my;i++){
            for(j=nx;j<ny;j++){
                di=M-1-(my-i);
                dj=M-1-(ny-j);
                A[i*order+j]=D[di*(M-1)+dj];
            }
        }
        for (i=nx;i<ny;i++){
            for(j=mx;j<my;j++){
                di=M-1-(ny-i);
                dj=M-1-(my-j);
                A[i*order+j]=D[di*(M-1)+dj];
            }
        }
    }
    int cx,cy;
    int ci,cj;

    for(n=0;n<N-1;n++){
        for(m=0;m<M-1;m++){
            C[m*(M-1)+m]=2.0*(Rx+Ry)+((m-M2+1)*hy)*((m-M2+1)*hy)+((n-N2+1)*hx)*((n-N2+1)*hx);
        }
        for(m=0;m<M-2;m++){
            C[m*(M-1)+m+1]=-Rx;
            C[(m+1)*(M-1)+m]=-Rx;
        }
        cx = (M-1)*n;
        cy = (M-1)*n+M-1;
        for (i=cx;i<cy;i++){
            for(j=cx;j<cy;j++)
            {
                ci=M-1-(cy-i);
                cj=M-1-(cy-j);
                A[i*order+j]=C[ci*(M-1)+cj];
            }
        }
    }

  int info = LAPACKE_dsyev(matrix_layout,jobz,uplo,order,A,lda,w);
  if(info==0)
    {
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
    mkl_free(D);
    mkl_free(C);

    return 0;
}
