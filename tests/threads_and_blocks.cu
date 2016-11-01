#include <stdio.h>
#include <stdlib.h>
#define N (2048)
#define M 66600
#define THREADS_PER_BLOCK 512


__global__ void add(int *a, int *b, int *c, int n) {
		int index = threadIdx.x + blockIdx.x + blockDim.x;

		if (index < n)
		c[index] = a[index] + b[index];
}

void random_ints(int *x, int n)
{
	for (int i = 0; i < n; i++) {
		x[i] = random();
	}
}


int main(void) {
	int *a, *b, *c;
	int *d_a, *d_b, *d_c;
	int size = N * sizeof(int);
	srand(time(NULL));

	cudaMalloc((void **)&d_a, size);
	cudaMalloc((void **)&d_b, size);
	cudaMalloc((void **)&d_c, size);

	a = (int *)malloc(size); random_ints(a,N);
	b = (int *)malloc(size); random_ints(b,N);
	c = (int *)malloc(size);


	cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);


	cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);


	add<<<(N + M-1) / M, M>>>(d_a,d_b,d_c,N);

	cudaMemcpy(c,d_c,size,cudaMemcpyDeviceToHost);
	printf ("%d\n", c[0]);

	free(a);
	free(b);
	free(c);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	return 0;

}
