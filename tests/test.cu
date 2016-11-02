#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

__device__ void dobra(int num)
{
	printf ("thread x2", num*2);
}

__global__ void sumAll(int *total)
{
	__shared__ int shared;
	shared+=threadIdx.x;
	dobra(threadIdx.x);

	__syncthreads();

	*total = shared;
}



int main(void) {

	int total =0;
	int *totalPtr;

	cudaMalloc((void **)&totalPtr, sizeof(int));
	cudaMemcpy(totalPtr, &total, sizeof(int), cudaMemcpyHostToDevice);

	sumAll<<<1,1024>>>(totalPtr);

	cudaMemcpy(&total, totalPtr, sizeof(int), cudaMemcpyDeviceToHost);

	printf("\nTotal: %i", total);

	cudaFree(totalPtr);

	return 0;
}
