#include <stdio.h>
#include <stdlib.h>

#include "MatUtil.h"

//mpicc -std=c99 -o APSPtest.c MathUtil.h MathUtil.c
//qsub -pe mpich 4 mysge.sh

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	if(argc != 2)
	{
		printf(" Usage: test {N}\n");
		exit(-1);
	}
	
	printf("Generating matrix...");
	//generate a random matrix.
	size_t N = atoi(argv[1]);
	int *mat = (int*)malloc(sizeof(int)*N*N);
	
	GenMatrix(mat, N);
	
	printf("Computing reference result...");
	//compute the reference result.
	int *ref = (int*)malloc(sizeof(int)*N*N);
	memcpy(ref, mat, sizeof(int)*N*N);
	ST_APSP(ref, N);
	
	printf("Running parallel algorithm...");
	
	//compute your results
	int *result = (int*)malloc(sizeof(int)*N*N);
	memcpy(result, mat, sizeof(int)*N*N);
	//replace by parallel algorithm
	MT_APSP(result, N);
	MPI_Finalize();

	printf("Parallel algorithm done.");
	
	//compare your result with reference result
	if(CmpArray(result, ref, N*N))
		printf("Your result is correct.\n");
	else
		printf("Your result is wrong.\n");

	
}
