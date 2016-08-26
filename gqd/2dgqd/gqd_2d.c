/********This is a code for two dimension g q d  eigen-problem************/
#include <stdio.h>
//lapacke headers
#include "mkl.h"

double POT(double x, double y)
{
    return -x*x-y*y;
}

int main(){
    int matrix_layout = LAPACK_COL_MAJOR;   //Careful!!
    char jobz = 'V';
    char uplo='U';
    int M,M2,N,N2;
    double xmax,xmin,ymax,ymin;//(range of space)
    double hx,hy,Rx,Ry;
    double out;
    M = 50;             // x direction divide into M slice
    N = 50;             // y direction divide into N slice
    M2 = M/2;
    N2 = N/2;
    xmax = 4.0;        // initial range of x, y;
    xmin = -4.0;
    ymax = 4.0;
    ymin = -4.0;
    
    hx = (xmax-xmin)/((double) N);      // step size of x direction;
    hy = (ymax-ymin)/((double) M);      // step size of y direction;
    Rx = 1.0/(hx);                   
    Ry = 1.0/(hy);
    int lda = 2*(M-1)*(N-1);           
    int order = 2*(M-1)*(N-1);
    MKL_Complex16 *A,*D,*C1,*C2;            // All Matrix labeled in note.
    double *eiv;
    MKL_Complex16 U[4],P[4],Q[4];       //     
    MKL_Complex16 M1[4],N1[4];
    FILE *s,*v;
    s=fopen("eigva.dat","w");
    v=fopen("eigve.dat","w");

    A = (MKL_Complex16 *)mkl_malloc( order*order*sizeof( MKL_Complex16 ), 64);
    eiv = (double *)mkl_malloc( order*sizeof( double ), 64);
    D = (MKL_Complex16 *)mkl_malloc( 2*(M-1)*2*(M-1)*sizeof( MKL_Complex16 ), 64);
    C1 = (MKL_Complex16 *)mkl_malloc( 2*(M-1)*2*(M-1)*sizeof( MKL_Complex16 ), 64);
    C2 = (MKL_Complex16 *)mkl_malloc( 2*(M-1)*2*(M-1)*sizeof( MKL_Complex16 ), 64);

    int i,j,m,n,mx,my,nx,ny,cx,cy;
    int di,dj,ci,cj;

/************************************************************/
/*                   initial M1,N1,P,Q                      */
/*             Attential matrix_layout   COL                */
    for(i=0;i<4;i++)
    {
        M1[i].real=M1[i].imag= N1[i].real=N1[i].imag= P[i].real=P[i].imag= Q[i].real=Q[i].imag=0.;
        U[i].real=U[i].imag=0.;
    }
    M1[1].imag=M1[2].imag=-Rx/2.;
    N1[1].imag=N1[2].imag=Rx/2.;
    P[1].real=Q[2].real=Ry/2.;
    P[2].real=Q[1].real=-Ry/2.;
/*                                                          */
/************************************************************/
	printf("initialized matrix M,N,P,Q......\n");
/************************************************************/
/*                   initial C1,C2                          */
for(m=0;m<M-1;m++){
        cx = 2*m;
        cy = 2*m+2;
        for (i=cx;i<cy;i++){
            for(j=cx;j<cy;j++)
            {
                ci=2-(cy-i);
                cj=2-(cy-j);
                C1[i*2*(M-1)+j]=P[ci*2+cj];
                C2[i*2*(M-1)+j]=Q[ci*2+cj];
            }
        }
    }
/*                                                          */
/************************************************************/

/************************************************************/
/*                 initial D and then A                     */
/*                A is the target Matrix                    */
/*          off-diagonalization terms of D Matrix           */
for(m=0;m<M-2;m++)
{
     mx = 2*m;
     my = 2*m+2;
     nx = 2*m+2;
     ny = 2*m+2+2;
    for (i=mx;i<my;i++){
        for(j=nx;j<ny;j++){
            di=2-(my-i);
            dj=2-(ny-j);
            D[i*2*(M-1)+j]=N1[di*2+dj];
        }
    }
    for (i=nx;i<ny;i++){
        for(j=mx;j<my;j++){
            di=2-(ny-i);
            dj=2-(my-j);
            D[i*2*(M-1)+j]=M1[di*2+dj];
        }
    }
}
printf("initialized matrix C1,C2,and off-diag D ......\n");
for(n=0;n<N-1;n++)
{
    /******** This circle is to initialize D *********/
    for(m=0;m<M-1;m++) // m, n stands for space coor.
    {
        /* initialize u Matrix */
        U[3].real=U[0].real=POT(xmin+(m+1)*hx,ymin+(n+1)*hy);
        /*  D diagonalization terms -- U matirx*/
        cx = 2*m;
        cy = 2*m+2;
        for (i=cx;i<cy;i++){
            for(j=cx;j<cy;j++)
            {
                ci=2-(cy-i);
                cj=2-(cy-j);
                D[i*2*(M-1)+j]=U[ci*2+cj];
            }
        }
    }
    /********   After ini-D it's to initialize A  use D and C1 and C2 *********/
    /*  A diagonalization terms -- D matirx*/
    cx = 2*(M-1)*n;
    cy = 2*(M-1)*n+2*(M-1);
    for (i=cx;i<cy;i++){
    for(j=cx;j<cy;j++)
    {
        ci=2*(M-1)-(cy-i);
        cj=2*(M-1)-(cy-j);
        A[i*order+j]=D[ci*(2*(M-1))+cj];
    }
    }
    /*  A off-diagonalization terms -- C1,C2 Matrix */
    if(n<N-2)
    {
    mx = 2*(M-1)*n;
    my = 2*(M-1)*n+2*(M-1);
    nx = 2*(M-1)*n+2*(M-1);
    ny = 2*(M-1)*n+2*(M-1)+2*(M-1);
    for (i=mx;i<my;i++){
        for(j=nx;j<ny;j++){
            di=2*(M-1)-(my-i);
            dj=2*(M-1)-(ny-j);
            A[i*order+j]=C2[di*2*(M-1)+dj];
            }
        }
    for (i=nx;i<ny;i++){
        for(j=mx;j<my;j++){
            di=2*(M-1)-(ny-i);
            dj=2*(M-1)-(my-j);
            A[i*order+j]=C1[di*2*(M-1)+dj];
        }
    }
    }
}
printf("initialized the target matrix A ......\n");
/*                                                          */
/************************************************************/


 int info = LAPACKE_zheev(matrix_layout,jobz,uplo,order,A,lda,eiv);

if(info==0)
{
    for(i=0;i<order;i++){
    //  fprintf(s,"eigenvalue %d:\n",i);
        fprintf(s,"%.6g\n",eiv[i]);       //print eigenvalue.
    // fprintf(v,"right eigenvector: \n");
        for(j=0;j<order;j+=2)
            {
    
                out=(A[i*order+j].real)*(A[i*order+j].real)+(A[i*order+j].imag)*(A[i*order+j].imag)+
                (A[i*order+j+1].real)*(A[i*order+j+1].real)+(A[i*order+j+1].imag)*(A[i*order+j+1].imag);
                fprintf(v,"%.6g\t",out);  //print eigenvector.
            }
            fprintf(v,"\n");
        }
        printf("SUCCESS\n");
    }
    fclose(s);
    fclose(v);

    mkl_free(A);
    mkl_free(eiv);
    mkl_free(D);
    mkl_free(C1);
    mkl_free(C2);

    return 0;
}
