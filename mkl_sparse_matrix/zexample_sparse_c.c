/*******************************************************************************
* Copyright 2005-2016 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
!   Content: Example for Intel MKL Extended Eigensolvers (sparse format,
!            double complex precision)
!
!*******************************************************************************
!
!
! The following routines are used in the example:
!          ZGEMM  ZFEAST_HCSREV  ZFEAST_HCSRGV  FEASTINIT.
!
! Consider the matrix A
!                 |  10   1+2i  0    0    0    0    0    0    0    0    |
!                 |  1-2i  9   2+3i  0    0    0    0    0    0    0    |
!                 |  0    2-3i  8   3+4i  0    0    0    0    0    0    |
!                 |  0    0    3-4i  7   4+5i  0    0    0    0    0    |
!                 |  0    0    0    4-5i  6   5+6i  0    0    0    0    |
!    A    =       |  0    0    0    0    5-6i  5   6+7i  0    0    0    |,
!                 |  0    0    0    0    0    6-7i  4   7+8i  0    0    |
!                 |  0    0    0    0    0    0    7-8i  3   8+9i  0    |
!                 |  0    0    0    0    0    0    0    8-9i  2   9+10i |
!                 |  0    0    0    0    0    0    0    0    9-10i 1    |
!
! stored as sparse matrix.
! B is a unit matrix:
!                 |  1   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0  |
!                 |  0   0   0   0   0   0   0   0   0   1  |
!
!  In what follows the symbol ' represents a conjugate transpose operation.
!
!  The test performs the following operations:
!
!  Step 1. Calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!  Step 2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using
!          ZFEAST_HCSREV.
!
!  Step 3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!           are the expected eigenvalues  and E(i) are eigenvalues computed
!           by ZFEAST_HCSREV().
!
!  Step 4. The code computes the maximum absolute value of elements
!          of the matrix Y = X' *X - I, where X is the matrix of eigenvectors
!          computed by ZFEAST_HCSREV.
!          ZGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!  Step 5. The  code solves  the generalized eigenvalue problem Ax=eBx using
!          ZFEAST_HCSRGV.
!
!  Step 6. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!           are the expected eigenvalues  and E(i) are eigenvalues computed
!           by ZFEAST_HCSRGV().
!
!  Step 7. The code computes the maximum absolute value of the elements of
!          the matrix  Y = X' * X - I, where X is the matrix of eigenvectors
!          computed by ZFEAST_HCSRGV.
!          ZGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mkl_solvers_ee.h"

#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    char          UPLO = 'F'; /* Type of matrix: (F=full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. We use size of matrix N and 3 arrays to store matrix in CSR format */
    const MKL_INT N = 10;
    MKL_INT       rows[11] = { 1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 29 };
    MKL_INT       cols[28] = {          1,            2,
                                        1,            2,            3,
                                                      2,            3,            4,
                                                                    3,            4,            5,
                                                                                  4,            5,            6,
                                                                                                5,            6,            7,
                                                                                                              6,            7,            8,
                                                                                                                            7,            8,            9,
                                                                                                                                          8,            9,           10,
                                                                                                                                                        9,           10
                            };
    MKL_Complex16 val[28] = {{10.0,  0.0}, { 1.0,  2.0},
                             { 1.0, -2.0}, { 9.0,  0.0}, { 2.0,  3.0},
                                           { 2.0, -3.0}, { 8.0,  0.0}, { 3.0,  4.0},
                                                         { 3.0, -4.0}, { 7.0,  0.0}, { 4.0,  5.0},
                                                                       { 4.0, -5.0}, { 6.0,  0.0}, { 5.0,  6.0},
                                                                                     { 5.0, -6.0}, { 5.0,  0.0}, { 6.0,  7.0},
                                                                                                   { 6.0, -7.0}, { 4.0,  0.0}, { 7.0,  8.0},
                                                                                                                 { 7.0, -8.0}, { 3.0,  0.0}, { 8.0,  9.0},
                                                                                                                               { 8.0, -9.0}, { 2.0,  0.0}, { 9.0, 10.0},
                                                                                                                                             { 9.0,-10.0}, { 1.0,  0.0}
                            };
     

    /* Declaration of FEAST variables */
    MKL_INT       fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */
    double        Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */

    double        epsout;        /* Relative error of the trace */
    MKL_INT       loop;          /* Number of refinement loop */
    MKL_INT       L = 8;
    MKL_INT       M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT       M;             /* Total number of eigenvalues found in the interval */

    double        E[10];         /* Eigenvalues */
    MKL_Complex16 X[100];        /* Eigenvectors */
    double        res[10];       /* Residual */
    /* Declaration of local variables */
    MKL_INT       info;          /* Errors */
    double        Eig[10];       /* Eig - array for storing exact eigenvalues */
    double        R[10];         /* R = |E-Eig| */
    MKL_Complex16 Y[10][10];     /* Y=(X')*X-I */

    char          ZGEMMC = 'C';  /* Character for GEMM routine, conjugated transposed case */
    char          ZGEMMN = 'N';  /* Character for GEMM routine, non-transposed case */
    MKL_Complex16 one  = {1.0, 0.0};    /* alpha parameter for GEMM */
    MKL_Complex16 zero = {0.0, 0.0};    /* beta  parameter for GEMM */


    MKL_INT      i, j;
    double        trace, smax, eigabs;
	clock_t start,stop;
	start=clock();

    for (i=4; i<N; i++)
    {
        Eig[i] = 0.0;
    }

    printf("\n FEAST ZFEAST_HCSREV AND ZFEAST_HCSRGV EXAMPLE\n");
    /* Initialize matrix X */
    for (i=0; i<N*N; i++)
    {
        X[i] = zero;
    }

    printf("Sparse matrix size %i\n", (int)N);

    /* Search interval [Emin,Emax] */
    Emin = 2.0;
    Emax = 12.0;
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

    fpm[1] =  48; /* Extended Eigensolver routines print runtime status to the screen. */

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
        printf("Routine zfeast_hcsrev returns code of ERROR: %i", (int)info);
        return 1;
    }

 for (i=0;i<M;i++)
 {
     printf("%d %f : ",i,E[i]);
     for (j=0;j< N;j++)
     {
        
         printf("{%f,%f}\t",X[i*N+j].real,X[i*N+j].imag);
   
    }
     printf("\n");
 }
    
    return 0;
}
