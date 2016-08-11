/****   This a code to cal magnititude of complex vector  ****/
#include <stdio.h>
#include <stdlib.h>
#include <mkl_cblas.h>

#define N	10
void initVector(MKL_Complex16 * v)
{
	int i;
	for (i = 0; i < N; i++) {
		v[i].real = -i * 1.0f;
		v[i].imag = i * 1.0f;
	}
}

int	main(int argc, char *argv[])
{
	MKL_Complex16  *vector = (MKL_Complex16 *) malloc(sizeof(MKL_Complex16) * N);
	initVector(vector);

	double	ret3 = cblas_dzasum(N, vector, 1); //cal magnititude. return a double number.
		printf("Result of sasum: %lf\n", ret3);
}
