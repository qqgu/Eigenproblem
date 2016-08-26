#include <stdio.h>
#include "mkl.h"
int main()
{
	MKL_Complex16  a[4];
	a[1].real=1;
	a[1].imag=1;
	MKL_Complex16 A[4];
//	A[0].real=0;
	double x;
	A[1].imag=0.1;
	int i;
//	a[0]=a[1]=a[2]=a[3]=1;
	A[1]=a[1];
	x=a[1]*a[1];
	for(i=0;i<4;i++)
	{
	printf("%f\t",x);
	}
	printf("\n");
}
