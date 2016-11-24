#include <stdio.h>
#include <unistd.h>
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
#define FREQUENCY_MAX = 100
#define FREQUENCY_MIN = 0

struct bat {
    //tends towards 1
    double pulse_rate;
    //tends towards 0
    double loudness;
    double fitness;
    double frequency;
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];
};


__device__ void update_velocity(struct bat *bat, struct bat *best)
{
    for (int i = 0; i < DIMENSIONS; ++i) {
        bat->velocity[i]+= (bat->position[i] - best->position[i]) * bat->frequency;
        /* printf("Velocity: %E\n", bat->velocity[i]); */
        /* force_boundry_on_value(&bat->velocity[i]); */
    }
} 

/* double (*objective_function)(double[], int); */

__device__ double sphere (double *solution, int dimensions)
{
    double total = 0;

    for (int i = 0; i < dimensions; i++) {
        total+= solution[i] * solution[i];
    }

    return total;
}

__device__ void copy_bat(struct bat *from, struct bat *to)
{
    memcpy(to, from, sizeof(struct bat));
}

__device__ void  get_best(struct bat *bats, struct bat *best)
{
    double current_best_val; 
    int best_indice;

    current_best_val = bats[0].fitness;
    best_indice = 0;
    for (int i = 0; i < BATS_COUNT; i++) {
        if (bats[i].fitness < current_best_val) {
            current_best_val = bats[i].fitness;
            best_indice = i;
        }
    }
    copy_bat(&bats[best_indice], best);
}


__device__ void log_bat_stdout(struct bat *bat, int dimensions) 
{
    /* printf("Best BAT \n"); */
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        /* printf("Position Dimension [%d]: %f\n", i, bat->position[i]); */
        position_average+=bat->position[i];
    }
    /* position_average/=dimensions; */
    printf("Frequency: %E\n", bat->frequency);
    printf("Pulse-rate: %E\n", bat->pulse_rate);
    printf("Loudness: %E\n", bat->loudness);
    printf("Position Average: %f\n", position_average);
    printf("Fitness E: %E\n", bat->fitness);
    /* printf("Fitness F: %f\n", bat->fitness); */
}

/* __device__ void update_position(struct bat *bat) */
/* { */
/*     for (int i = 0; i < DIMENSIONS; ++i) { */

/*         /1* force_boundry_on_value(&bat->position[i]); *1/ */
/*     } */
/* } */


__global__ void run_bats(unsigned int seed, struct bat *bats, struct bat *candidates)
{
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    __shared__ struct bat *best;
    __shared__ int iteration;
    best = (struct bat *) malloc(sizeof(struct bat));


    bats[threadIdx.x].pulse_rate = 0.0;
    bats[threadIdx.x].frequency = 0.0;
    bats[threadIdx.x].fitness;
    bats[threadIdx.x].loudness = INITIAL_LOUDNESS;
    for (int j = 0; j < DIMENSIONS; j++) {
        bats[threadIdx.x].velocity[j] = 0;
        bats[threadIdx.x].position[j] = curand(&state);
    }

    bats[threadIdx.x].fitness = sphere(bats[threadIdx.x].position, DIMENSIONS);

    __syncthreads();

    get_best(bats, best);

    iteration = 0;
    while(iteration < ITERATIONS) {
            //frequency
            double beta = curand_uniform(&state);
            bats[threadIdx.x].frequency = beta;
            update_velocity(&bats[threadIdx.x], best);
            copy_bat(&bats[threadIdx.x], &candidates[threadIdx.x]);

            for (int i = 0; i < DIMENSIONS; ++i) {
                candidates[threadIdx.x].position[i] += candidates[threadIdx.x].velocity[i];
            }

            //local search
            if (curand_uniform(&state) < candidates[threadIdx.x].pulse_rate) {
                for (int i = 0; i < DIMENSIONS; i++ ) {
                    candidates[threadIdx.x].position[i] = best->position[i] +  curand_uniform(&state);
                }
            }


            bats[threadIdx.x].fitness = sphere(bats[threadIdx.x].position, DIMENSIONS);
            candidates[threadIdx.x].fitness = sphere(candidates[threadIdx.x].position, DIMENSIONS);

            if (candidates[threadIdx.x].fitness < bats[threadIdx.x].fitness) {
                printf("fitness");
                copy_bat(&candidates[threadIdx.x], &bats[threadIdx.x]);
                bats[threadIdx.x].pulse_rate = 1 - exp(-LAMBDA*iteration);
                bats[threadIdx.x].loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);

            }

            __syncthreads();
            iteration++;
            get_best(bats, best);
    }

    /* log_bat_stdout(best, DIMENSIONS); */
}


int main(void)
{
    struct bat *bats;
    struct bat *candidates;
    struct bat *best;
    cudaMalloc((void **)&bats, sizeof(struct bat) * BATS_COUNT);
    cudaMalloc((void **)&candidates, sizeof(struct bat) * BATS_COUNT);

    run_bats<<<1,BATS_COUNT>>>(time(NULL), bats, candidates);
    cudaDeviceSynchronize();

    cudaFree(bats);
    return 0;
}

