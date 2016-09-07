#include "MatUtil.h"
#include <mpi.h>

void GenMatrix(int *mat, const size_t N)
{
	for(int i = 0; i < N*N; i ++)
		mat[i] = rand()%32 - 1;
	for(int i = 0; i < N; i++)
		mat[i*N + i] = 0;

}

bool CmpArray(const int *l, const int *r, const size_t eleNum)
{
	for(int i = 0; i < eleNum; i ++)
		if(l[i] != r[i])
		{
			printf("ERROR: l[%d] = %d, r[%d] = %d\n", i, l[i], i, r[i]);
			return false;
		}
	return true;
}


/*
	Sequential (Single Thread) APSP on CPU.
*/
void ST_APSP(int *mat, const size_t N)
{
	for(int k = 0; k < N; k ++)
		for(int i = 0; i < N; i ++)
			for(int j = 0; j < N; j ++)
			{
				int i0 = i*N + j;
				int i1 = i*N + k;
				int i2 = k*N + j;
				if(mat[i1] != -1 && mat[i2] != -1)
                     { 
			          int sum =  (mat[i1] + mat[i2]);
                          if (mat[i0] == -1 || sum < mat[i0])
 					     mat[i0] = sum;
				}
			}
}

/*
Parallel (Multiple Threads) APSP on CPU.
*/
void MT_APSP(int *mat, const size_t N)
{
	MPI_Init(&argc, &argv);
	for (int k = 0; k < N; k++)
	{
		
		functionName(int *m);
		
	}
	MPI_Finalize();
}

void functionName(int *mat, const size_t N, const size_t K)
{

	int rank, numprocs;

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int lb = (N / numprocs)*rank;
	int ub = (N / numprocs)*(rank+1);
	
	for (int i = lb; i < ub; i++)
	{
		for (int j = 0; j < N; j++)
		{
			int i0 = i*N + j;
			int i1 = i*N + k;
			int i2 = k*N + j;
			if (mat[i1] != -1 && mat[i2] != -1)
			{
				int sum = (mat[i1] + mat[i2]);
				if (mat[i0] == -1 || sum < mat[i0])
					mat[i0] = sum;
			}
		}
	}
}


