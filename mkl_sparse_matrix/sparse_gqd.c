/*
!   This code is designed to solve gqd eigenvalue proble 2-dimension
!   using mkl extended eigen solves
!   sparse matrix is stored using CSR format rows[0] = 1
!   Autor: Qiangqiang Gu
!   Time : 2016 - 08 - 31
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include <time.h>
#include "mkl_solvers_ee.h"


double POT (double x,double y)
{
    return -x*x - y*y;
}

 
 int main()
 {
/*
!   Define some variables of specific problem.
!   range of x, and y, size of origin matrix, matrix block;
*/
    MKL_INT Nd = 50;                          // number of discrete
    MKL_INT Ncount = Nd-1;
    MKL_INT Nblc = 2*Ncount;                //matrix block; 
    const MKL_INT N =  2*Ncount*Ncount;     //origin matrix is N by N matirx.

    double xmin = -4.0;
    double xmax = 4.0;
    double ymin = -4.0;
    double ymax = 4.0;
    double hx,hy,Rx,Ry;                    //delta x = hx; delta y = hy, R = 1/2h;
    
    hx = (xmax-xmin)/((double) Nd);         // step size of x direction;
    hy = (ymax-ymin)/((double) Nd);         // step size of y direction;
    Rx = 1.0/(2.0*hx);                   
    Ry = 1.0/(2.0*hy);

/*
!   use CSR format to store matrix
!   rows, cols , and vals arrays is the three arrays that CSR needed
*/   

   
    MKL_INT rows[N+1];                              
    MKL_INT rval[Nblc];                              //use rval to init rows;
    //rows = (MKL_INT *)mkl_malloc( (N+1)*sizeof( MKL_INT ), 64);
    //rval = (MKL_INT *)mkl_malloc( (Nblc)*sizeof( MKL_INT ), 64);
    
    // some assistant variables; cycle etc.
    int i,j,m,flag,ic,flag2;
    int numr;
    double mod;
    FILE *s,*v;
    s=fopen("eigva.dat","w");
    v=fopen("eigve.dat","w");
/*
!   init rows array.
!   first elements of rows rows[0] = 1; is requaired by fast extended routines! should be careful.
!   rval is a assistant array to initialize rows.   
*/
	clock_t start,stop;
	start=clock();

    rows[0]=1;
    for(i=0;i<Nblc;i+=2)
    {
        if(i==0||i==Nblc-2)
        {
            rval[i+1]=rval[i]=3;
        } 
        else rval[i+1]=rval[i]=4;       
    }
    for (i=0;i<N;i++)
    {  
        if (i<Nblc) rows[i+1]=rows[i]+rval[i];
        else if(i<(Ncount-1)*Nblc) rows[i+1]=rows[i]+rval[i%(Nblc)]+1;
        else rows[i+1]=rows[i]+rval[i%(Nblc)]; 
    }

/*
!      init val matrix all the nonzero elements of origin matrix M    
!              rows[N]-1 is the number of nonzero elements    
!      init cols of all nonzero elements col index; during init val.     
*/
    MKL_INT cols[rows[N]-1];                          // nonzero elements col index
    MKL_Complex16 val[rows[N]-1];                 // nonzero elements value
    MKL_INT u[Nblc];                             // u array is store the order of every row nonzero element sequence.

   // cols = (MKL_INT *)mkl_malloc( (rows[N]-1)*sizeof( MKL_INT ), 64);
    //val = (MKL_Complex16 *)mkl_malloc( (rows[N]-1)*sizeof( MKL_Complex16 ), 64);
    //u = (MKL_INT *)mkl_malloc( (Nblc)*sizeof( MKL_INT ), 64);
 
 /**************  initialize val and cols part which is the most tricky part of this code.  ******************/

    u[0]=0;
    u[1]=0;
    for(i=2;i<Nblc;i++)  //for the first block u is {0,0,1,1...1,1}; and others is {1,1,2,2...2,2}
    {
        u[i]=1;
    }

    for (i=0;i<N;i++)
    {
        j = i/Nblc;
        ic = i%Nblc;
        flag2 = i%2;
        numr = rows[i+1]-rows[i];   // nonzero elements number of row i; 
        if(i < Nblc)    //the first block has only D and C1 matrix block
        {
            val[rows[i]+u[ic]-1].real = POT(xmin+(ic/2+1)*hx,ymin+(j+1)*hy);
            val[rows[i]+u[ic]-1].imag = 0;
            
            cols[rows[i]+u[ic]-1] = i+1;  // diagonal element i+1 for first index is 1 not 0

           for( m=u[ic]+1;m<numr-1;m++)
            {   
                    val[rows[i]+m-1] = (MKL_Complex16) {0.0,-1.0*Rx};
                    cols[rows[i]+m-1] = i+1+2+pow(-1,flag2);
            }

            val[rows[i+1]-1-1] = (MKL_Complex16) {pow(-1,flag2)*(-1)*Ry,0.0}; 
            cols[rows[i+1]-1-1] = i+1+Nblc+pow(-1,flag2);
        
            for(m=u[ic]-1;m>=0;m--)
            {
                val[rows[i]+m-1] = (MKL_Complex16) {0.0,Rx};
                cols[rows[i]+m-1] = i+1-(2-pow(-1,flag2));
            }
        }
        
        if(i>=Nblc&&i<(Ncount-1)*Nblc)     //middle blocks have C2 D and C1 blocks.
        {
            val[rows[i]+u[ic]+1-1].real = POT(xmin+(ic/2+1)*hx,ymin+(j+1)*hy);
            val[rows[i]+u[ic]+1-1].imag = 0;
            cols[rows[i]+u[ic]+1-1] = i+1;
    
            for( m=u[ic]+1+1;m<numr-1;m++)
            {   
                val[rows[i]+m-1] = (MKL_Complex16) {0,-1.0*Rx};
                cols[rows[i]+m-1] = i+1+2+pow(-1,flag2);
            }
            
            val[rows[i+1]-1-1] = (MKL_Complex16) {pow(-1,flag2)*(-1)*Ry,0}; 
            cols[rows[i+1]-1-1] = i+1+Nblc+pow(-1,flag2);
            

            for(m=u[ic];m>0;m--)    
            {
                val[rows[i]+m-1] = (MKL_Complex16) {0,Rx};
                cols[rows[i]+m-1] = i+1-(2-pow(-1,flag2));
            }

            val[rows[i]-1] = (MKL_Complex16) {pow(-1,flag2)*Ry,0}; 
            cols[rows[i]-1] = i+1-(Nblc-pow(-1,flag2));

  
        }
        if(i>=(Ncount-1)*Nblc)  //The last block has C2 and D matrix block
        {
            val[rows[i]+u[ic]+1-1].real = POT(xmin+(ic/2+1)*hx,ymin+(j+1)*hy);
            val[rows[i]+u[ic]+1-1].imag = 0;
            cols[rows[i]+u[ic]+1-1] = i+1;
            
            for( m=u[ic]+1+1;m<numr;m++)
            {   
                val[rows[i]+m-1] = (MKL_Complex16) {0,-1.0*Rx};
                cols[rows[i]+m-1] = i+1+2+pow(-1,flag2);
            }


            for(m=u[ic];m>0;m--)
            {
                val[rows[i]+m-1] = (MKL_Complex16) {0,Rx};
                cols[rows[i]+m-1] = i+1-(2-pow(-1,flag2));
            }

            val[rows[i]-1] = (MKL_Complex16) {pow(-1,flag2)*Ry,0}; 
            cols[rows[i]-1] = i+1-(Nblc-pow(-1,flag2));

        }


    }

/* Declaration of FEAST routine  variables */
    char        UPLO = 'F';     /* Type of matrix: (F=full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. We use size of matrix N and 3 arrays to store matrix in CSR format */

    MKL_INT     fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */
    double      Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */
    double      epsout;        /* Relative error of the trace */
    MKL_INT     loop;           /* Number of refinement loop */
    MKL_INT     L = 500;       /* guess of number of Eigenvalues between Emin and Emax.*/
    MKL_INT     M0 = L;         /* Initial guess for subspace dimension to be used */
    MKL_INT     M  = L;              /* Total number of eigenvalues found in the interval */

    double        E[L];          /* Eigenvalues */
    MKL_Complex16 X[L*N];        /* Eigenvectors */
    double        res[L];        /* Residual */

    MKL_INT     info;
    MKL_Complex16 zero = {0.0, 0.0};
    
    /* Initialize matrix X */
    for (i=0; i<L*N; i++)
    {
        X[i] = zero;
    }
    Emax = -2.0;
    Emin = -30.0;

    loop = 0;
    info = 0;
    epsout = 0.0;


    /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
    feastinit(fpm);     /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
    
    /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
    zfeast_hcsrev(
        &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
        &N,      /* IN: Size of the problem */
        val,     /* IN: CSR matrix A, values of non-zero elements */
        rows,    /* IN: CSR matrix A, index of the first non-zero element in row */
        cols,    /* IN: CSR matrix A, columns indices for each non-zero element */
        fpm,     /* IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
        &epsout, /* OUT: Relative error of on the trace */
        &loop,   /* OUT: Contains the number of refinement loop executed */
        &Emin,   /* IN: Lower bound of search interval */
        &Emax,   /* IN: Upper bound of search interval */
        &M0,     /* IN: The initial guess for subspace dimension to be used. */
        E,       /* OUT: The first M entries of Eigenvalues */
        X,       /* IN/OUT: The first M entries of Eigenvectors */
        &M,      /* OUT: The total number of eigenvalues found in the interval */
        res,     /* OUT: The first M components contain the relative residual vector */
        &info    /* OUT: Error code */
        );

    if ( info != 0 )
    {
        printf("Routine zfeast_hcsrev returns code of ERROR: %i\n", (int)info);
        return 1;
    }
    else printf ("SUCCESS!\n");
        
/*
！   output control
!   FEAST routine is col main order be careful to this！
!   just print out the moudle of the wavefuncion
!   |u1| + |u2| is |Psi|
!   an eigenvalue print out in one row to wave outfile.
!   each eigenvalue corresponding to a row in eig outfile 
*/


for(i=0;i < M;i++)
{
    fprintf(s,"%.6f\n",E[i]);       //print out eigenvalue to eigva.dat
    for (j = 0;j < N; j+=2)
    {
        mod=(X[i*N+j].real)*(X[i*N+j].real)+(X[i*N+j].imag)*(X[i*N+j].imag)+
           (X[i*N+j+1].real)*(X[i*N+j+1].real)+(X[i*N+j+1].imag)*(X[i*N+j+1].imag);
        fprintf(v,"%.6f\t",mod);  //print out eigenvectors to eigve.dat
    }
   fprintf(v,"\n");

}

stop=clock();
printf("Elapsed time = %.5f seconds\n",((double)(stop - start))/CLOCKS_PER_SEC);
fclose(s);
fclose(v);


//    mkl_free(rows);
 //   mkl_free(rval);
  //  mkl_free(cols);
   // mkl_free(val);
    //mkl_free(u);

return 0;
 }