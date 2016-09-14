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
	
	const int root = 0;
	int tmp_k;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == root)
	{
		n = N;
		printMat(mat, n);
	}
			
	//broadcast N
	MPI_Bcast (&n, 1, MPI_INT, root, MPI_COMM_WORLD);
	printf("rank = %d n = %d\n", rank, n);
	int *tmp_buf;
	tmp_buf = malloc(n*n*sizeof(int));
	
	
	for (int k = 0; k < N; k++)
	{
		if(rank == root)
		{
			tmp_k = k;
			for(int i = 0; i < n*n; i++)
				tmp_buf[i] = mat[i];
		}
		
		//printf("Rank %d broadcasting current k %d...\n", rank, k);
		//MPI_Bcast (&tmp_k, 1, MPI_INT, root, MPI_COMM_WORLD);
		
		printf("Rank %d broadcasting current matrix...\n", rank);
		MPI_Bcast (tmp_buf, n*n, MPI_INT, root, MPI_COMM_WORLD);
		
		
		
		printf("Rank %d calculating current matrix...\n", rank);
		functionName(tmp_buf, n, k );

		printf("Rank %d collecting data with reduce...\n", rank);
		MPI_Reduce(tmp_buf, mat, n*n, MPI_INT, MPI_MIN, root, MPI_COMM_WORLD);
		if(rank == 0)
			printMat(mat, n);
	    
	}
	
}

void functionName(int *mat, const int N, const int K)
{
	
	int rank, numprocs;
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int lb = (N / numprocs)*rank;
	int ub = (N / numprocs)*(rank+1);
	printf("rank = %d K = %d N = %d lb = %d ub = %d \n", rank, K, N, lb, ub);
	
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

void printMat(int *mat, const int N)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("salut rank = %d!\n", rank);
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			printf(" %d", mat[i*N + j]);
		}
		printf("\n");
	}
	printf("Au revoir!\n");	
}


