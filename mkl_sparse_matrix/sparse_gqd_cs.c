#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include "mkl_solvers_ee.h"

double POT (double x,double y)
{
    return -x*x - y*y;
}

 
 int main()
 {
    char  UPLO = 'F'; /* Type of matrix: (F=full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. We use size of matrix N and 3 arrays to store matrix in CSR format */
   
    MKL_INT Nd = 6; // number of discrete.
    MKL_INT Ncount = Nd-1;
    MKL_INT Nblc = 2*Ncount;    //matrix block; 
    const MKL_INT N =  2*Ncount*Ncount; //origin matrix is N by N matirx.

    double xmin = -4.0;
    double xmax = 4.0;
    double ymin = -4.0;
    double ymax = 4.0;
    double hx,hy,Rx,Ry;     //delta x = hx; delta y = hy, R = 1/h;
    
    hx = (xmax-xmin)/((double) Nd);      // step size of x direction;
    hy = (ymax-ymin)/((double) Nd);      // step size of y direction;
    Rx = 1.0/(2.0*hx);                   
    Ry = 1.0/(2.0*hy);
    
    int *rows;
    int *rval; //use rval to init rows;
    rows = (MKL_INT *)mkl_malloc( (N+1)*sizeof( MKL_INT ), 64);
    rval = (MKL_INT *)mkl_malloc( (Nblc)*sizeof( MKL_INT ), 64);
    
    int i,j,m,flag,ic,flag2;
    int numr;
/**                     init rows                **/
    for(i=0;i<Nblc;i+=2)
    {
        if(i==0||i==Nblc-2)
        {
            rval[i+1]=rval[i]=3;
        } 
        else rval[i+1]=rval[i]=4;       
    }

    rows[0]=0;

    for (i=0;i<N+1;i++)
    {  
        if (i<Nblc) rows[i+1]=rows[i]+rval[i];
        else if(i<(Ncount-1)*Nblc) rows[i+1]=rows[i]+rval[i%(Nblc)]+1;
        else rows[i+1]=rows[i]+rval[i%(Nblc)]; 
    }
/**       init val matrix all the nonzero elements of origin matrix M    **/
/*              rows[N] is the number of nonzero elements               */
    int *cols;  // nonzero elements col index
    MKL_Complex16 *val; // nonzero elements value
    int *u;
    cols = (MKL_INT *)mkl_malloc( (rows[N])*sizeof( MKL_INT ), 64);
    val = (MKL_Complex16 *)mkl_malloc( (rows[N])*sizeof( MKL_Complex16 ), 64);
    u = (MKL_INT *)mkl_malloc( (Nblc)*sizeof( MKL_INT ), 64);
 // envalue diagonal term u for it is change for diff i,j;
 
    u[0]=0;
    u[1]=0;
    for(i=2;i<Nblc;i++)
    {
        u[i]=1;
    }

    for (i=0;i<N;i++)
    {
        j = i/Nblc;
        ic = i%Nblc;
        flag2 = i%2;
        numr = rows[i+1]-rows[i];   // nonzero elements number of row i; 
        if(i < Nblc) 
        {
            val[rows[i]+u[ic]].real = POT(xmin+(ic/2+1)*hx,ymin+(j+1)*hy);
            val[rows[i]+u[ic]].imag = 0;
            
            cols[rows[i]+u[ic]] = i;

           for( m=u[ic]+1;m<numr-1;m++)
            {   
                    val[rows[i]+m] = (MKL_Complex16) {0.0,-1.0*Rx};
                    cols[rows[i]+m] = i+2+pow(-1,flag2);
            }

            val[rows[i+1]-1] = (MKL_Complex16) {pow(-1,flag2)*(-1)*Ry,0.0}; 
            cols[rows[i+1]-1] = i+Nblc+pow(-1,flag2);
        
            for(m=u[ic]-1;m>=0;m--)
            {
                val[rows[i]+m] = (MKL_Complex16) {0.0,Rx};
                cols[rows[i]+m] = i-(2-pow(-1,flag2));
            }
        }
        
        if(i>=Nblc&&i<(Ncount-1)*Nblc)
        {
            val[rows[i]+u[ic]+1].real = POT(xmin+(ic/2+1)*hx,ymin+(j+1)*hy);
            val[rows[i]+u[ic]+1].imag = 0;
            cols[rows[i]+u[ic]+1] = i;
    
            for( m=u[ic]+1+1;m<numr-1;m++)
            {   
                val[rows[i]+m] = (MKL_Complex16) {0,-1.0*Rx};
                cols[rows[i]+m] = i+2+pow(-1,flag2);
            }
            
            val[rows[i+1]-1] = (MKL_Complex16) {pow(-1,flag2)*(-1)*Ry,0}; 
            cols[rows[i+1]-1] = i+Nblc+pow(-1,flag2);
            

            for(m=u[ic];m>0;m--)
            {
                val[rows[i]+m] = (MKL_Complex16) {0,Rx};
                cols[rows[i]+m] = i-(2-pow(-1,flag2));
            }

            val[rows[i]] = (MKL_Complex16) {pow(-1,flag2)*Ry,0}; 
            cols[rows[i]] = i-(Nblc-pow(-1,flag2));

  
        }
        if(i>=(Ncount-1)*Nblc)
        {
            val[rows[i]+u[ic]+1].real = POT(xmin+(ic/2+1)*hx,ymin+(j+1)*hy);
            val[rows[i]+u[ic]+1].imag = 0;
            cols[rows[i]+u[ic]+1] = i;
            
            for( m=u[ic]+1+1;m<numr;m++)
            {   
                val[rows[i]+m] = (MKL_Complex16) {0,-1.0*Rx};
                cols[rows[i]+m] = i+2+pow(-1,flag2);
            }


            for(m=u[ic];m>0;m--)
            {
                val[rows[i]+m] = (MKL_Complex16) {0,Rx};
                cols[rows[i]+m] = i-(2-pow(-1,flag2));
            }

            val[rows[i]] = (MKL_Complex16) {pow(-1,flag2)*Ry,0}; 
            cols[rows[i]] = i-(Nblc-pow(-1,flag2));

        }


    }







 for(i=0;i<N;i++)
    {
       for (j=rows[i];j<=rows[i+1]-1;j++)
       {
            printf ("{%d: %f,%f}\n",cols[j],val[j].real,val[j].imag);
       }
  //printf ("%d\n",i/Nblc+1);
    printf("\n");

    }
 //   for (i=0;i<Nblc;i++)
    {
   //     printf("%d\t",u[i]);
    }
    mkl_free(rows);
    mkl_free(rval);
    mkl_free(cols);
    mkl_free(val);
    mkl_free(u);

return 0;
 }