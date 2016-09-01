#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mkl_solvers_ee.h"

#define max(a, b) (a) < (b) ? (b): (a)
//[1    0    3+i     0]
//[0    0    3-i       0]
//[3-i  3+i    0       0]
//[0    0    0       4]
int main()
{
    char          UPLO = 'F'; /* Type of matrix: (F=full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. We use size of matrix N and 3 arrays to store matrix in CSR format */
    const MKL_INT N = 4;
    MKL_INT       rows[5] = { 1, 3, 4, 6, 7 };
    MKL_INT       cols[6] = { 1, 3, 3, 1, 2, 4};
    MKL_Complex16 val[6] =  {{1,0}, {3,1}, {3,-1}, {3,-1},{3,1}, {4,0}};

    /* Declaration of FEAST variables */
    MKL_INT       fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */
    double        Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */

    double        epsout;        /* Relative error of the trace */
    MKL_INT       loop;          /* Number of refinement loop */
    MKL_INT       L = 4;
    MKL_INT       M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT       M;             /* Total number of eigenvalues found in the interval */

    double        E[4];         /* Eigenvalues */
    MKL_Complex16 X[N*N];        /* Eigenvectors */
    double        res[4];       /* Residual */
    /* Declaration of local variables */
    MKL_INT       info;          /* Errors */

    MKL_INT      i, j;
    double        trace, smax, eigabs;
	clock_t start,stop;
	start=clock();
    MKL_Complex16 zero={0,0};

    /* Initialize matrix X */
    for (i=0; i<N*N; i++)
    {
        X[i] = zero;
    }

    printf("Sparse matrix size %i\n", (int)N);

    /* Search interval [Emin,Emax] */
    Emin = -5.0;
    Emax = 5.0;
    printf("Search interval [ %.15e, %.15e  ]  \n", Emin, Emax);

    M0   = L;
    M    = L;
    loop = 0;
    info = 0;
    epsout = 0.0;

    /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
    feastinit(
        fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
        );

 //   fpm[0] =  0; /* Extended Eigensolver routines print runtime status to the screen. */

    /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
    printf(" Testing zfeast_hcsrev\n");
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
    printf("FEAST OUTPUT INFO %d \n",info);
    if ( info != 0 )
    {
        printf("Routine zfeast_hcsrev returns code of ERROR: %i\n", (int)info);
        return 1;
    }

 for (i=0;i<M;i++)
 {
  //  printf("%f : ",E[i]);
     for (j=0;j<N;j++)
     {
        double mode = X[i*N+j].real*X[i*N+j].real+X[i*N+j].imag*X[i*N+j].imag;
         printf("%f\t",mode);

    }
     printf("\n");
 }

    
    return 0;
}
