/********This is a code for radial graphene quantum dots eigen-problem************/
#include <stdio.h>
//lapacke headers
#include "mkl.h"

double U(double r)
{
    return -r*r;
}

int main(){
    int matrix_layout = LAPACK_COL_MAJOR;
    char jobz = 'V';
    char uplo='U';
    int M,M2,N;
    double rmax,rmin;//(range of space)
    double h,R,mqn=1.5;     //megnet quantum number.
    M = 1200;
    N = 2;
    M2 = M/2;
    rmax = 12.0;
    rmin = 0.0;
    
    h = (rmax-rmin)/((double) M);
    R = 1.0/(2.0*h);
    int lda = (M-1)*(N);
    int order = (M-1)*N;
    double *A,*w,*Dp,*Dn,*C;
    FILE *s,*v;
    s=fopen("eigva.dat","w");
    v=fopen("eigve.dat","w");

    A = (double *)mkl_malloc( order*order*sizeof( double ), 64);
    w = (double *)mkl_malloc( order*sizeof( double ), 64);
    Dp = (double *)mkl_malloc( N*N*sizeof( double ), 64);
    Dn = (double *)mkl_malloc( N*N*sizeof( double ), 64);
    C = (double *)mkl_malloc( N*N*sizeof( double ), 64);
    
    int i,j,m,n,mx,my,nx,ny;
    int di,dj;
    for (m=0;m<N-1;m++){
    
        Dp[m+1+m*N]=R;
        Dn[m+1+m*N]=-R; //0 row 1 col

        Dp[m+(m+1)*N]=-R;
        Dn[m+(m+1)*N]=R;  //1 row 0 col 
    


    }
    for(m=0;m<M-2;m++){
        mx = (N)*m;
        my = (N)*m+N;
        nx = (N)*m+N;
        ny = (N)*m+N+N;
        for (i=mx;i<my;i++){
            for(j=nx;j<ny;j++){
                di=N-(my-i);
                dj=N-(ny-j);
                A[i*order+j]=Dp[di*N+dj];
            }
        }
        for (i=nx;i<ny;i++){
            for(j=mx;j<my;j++){
                di=N-(ny-i);
                dj=N-(my-j);
                A[i*order+j]=Dn[di*N+dj];
            }
        }
    }
    int cx,cy;
    int ci,cj;

    for(m=0;m<M-1;m++){
        for(n=0;n<N;n++){
            C[n*N+n]=U((m+1)*h);
        }
        for(n=0;n<N-1;n++){
            C[n*N+n+1]=mqn/((m+1)*h);
            C[(n+1)*(N)+n]=mqn/((m+1)*h);
        }
        cx = N*m;
        cy = N*m+N;
        for (i=cx;i<cy;i++){
            for(j=cx;j<cy;j++)
            {
                ci=N-(cy-i);
                cj=N-(cy-j);
                A[i*order+j]=C[ci*N+cj];
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
    mkl_free(Dn);
    mkl_free(Dp);
    mkl_free(C);

    return 0;
}