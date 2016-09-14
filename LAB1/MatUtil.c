#include "MatUtil.h"

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
	
	
	int rank;
	int n;
	int *buf;//int buf[N];
	const int root = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == root)
	{
		buf = mat;
		n = N;
	}
		
	//broadcast N
	MPI_Bcast (n, 1, MPI_INT, root, MPI_COMM_WORLD);
	

	for (int k = 0; k < N; k++)
	{
		//broadcast k
		MPI_Bcast (k, 1, MPI_INT, root, MPI_COMM_WORLD);
		
		if(rank == root)
			buf = mat;
		//broadcast matrix
		MPI_Bcast (&buf, n, MPI_INT, root, MPI_COMM_WORLD);
		
		//do calculations for this k
		functionName(buf, n, k );
		
		//collect matrices with reduce
		MPI_Reduce(&buf, &mat, n, MPI_INT, MPI_MIN, root, MPI_COMM_WORLD);
	    
	}
	
}

void functionName(int *mat, const int N, const int K)
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
			int i1 = i*N + K;
			int i2 = K*N + j;
			if (mat[i1] != -1 && mat[i2] != -1)
			{
				int sum = (mat[i1] + mat[i2]);
				if (mat[i0] == -1 || sum < mat[i0])
					mat[i0] = sum;
			}
		}
	}
}


