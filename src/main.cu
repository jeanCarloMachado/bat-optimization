#include <stdio.h>
#include <unistd.h>
extern "C" {
#include "bat.h"
}

#include <curand.h>
#include <curand_kernel.h>


#define ITERATIONS 500
#define BATS_COUNT 40
#define INITIAL_LOUDNESS 1.0
#define DIMENSIONS 100

//probability of accepting bad results
#define ALFA 0.5
//affects local search
#define LAMBDA 0.1

#define BETA_MAX 1.0
#define BETA_MIN -1.0

#define BOUNDRY_MAX = 100
#define BOUNDRY_MIN = 0


__global__ void run_bats(unsigned int seed, struct bat *bats)
{
    bats[threadIdx.x].loudness = INITIAL_LOUDNESS;

    bats[threadIdx.x].velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
    bats[threadIdx.x].position = (double *) malloc(sizeof(double) * DIMENSIONS);

    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);


    for (int j = 0; j < DIMENSIONS; j++) {
        bats[threadIdx.x].velocity[j] = 0;
        bats[threadIdx.x].position[j] = curand(&state);
    }
}

int main(void) {


    struct bat *bats;
    cudaMalloc((void **)&bats, sizeof(struct bat) * BATS_COUNT);
    run_bats<<<1,BATS_COUNT>>>(time(NULL), bats);
    cudaDeviceSynchronize();

    cudaFree(bats);
    return 0;
}

