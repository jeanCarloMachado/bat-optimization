#include <stdio.h>
#include <unistd.h>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>
extern "C" {
#include "bat.h"
}

#define LAMBDA 0.1
#define ALFA 0.5
#define BETA_MAX 1.0
#define BETA_MIN 0.0
#define INITIAL_LOUDNESS 1.0
#define DIMENSIONS 100

int iterations = 10000;
int bats_count = 768;
__device__ int dbats_count;
__device__ int diterations;
const int EVALUTAION_FUNCTION = ACKLEY;

#define CUDA_CALL(cuda_function, ...)  { \
    cudaError_t status = cuda_function(__VA_ARGS__); \
    cudaEnsureSuccess(status, #cuda_function, false, __FILE__, __LINE__); \
}

bool cudaEnsureSuccess(cudaError_t status, const char* status_context_description,
        bool die_on_error, const char* filename, unsigned line_number) {
    if (status_context_description == NULL)
        status_context_description = "";
    if (status == cudaSuccess) {
        return true;
    }
    const char* errorString = cudaGetErrorString(status);
    fprintf(stderr, "CUDA Error: ");
    if (status_context_description != NULL) {
        fprintf(stderr, "%s\n", status_context_description);
    }
    if (errorString != NULL) {
        fprintf(stderr,"%s\n", errorString);
    } else {
        fprintf(stderr, "(Unknown CUDA status code %i", status);
    }

    fprintf(stderr, "Filename: %s, Line: %i\n", filename, line_number);

    if(die_on_error) {
        exit(EXIT_FAILURE);
    }
    return false;
}

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


__device__ int BOUNDRY_MAX;
__device__ int BOUNDRY_MIN;
__device__ int FREQUENCY_MIN;
__device__ int FREQUENCY_MAX;

__device__ double (*objective_function)(double[], int);
__device__ double griewank(double solution[], int dimensions)
{
    double total = 0;

    double top1=0;
    double top2=1;

    for(int i=0;i<dimensions;i++)
    {
        top1=top1+pow((solution[i]),(double)2);
        top2=top2*cos((((solution[i])/sqrt((double)(i+1)))*M_PI)/180);
    }
    total=(1/(double)4000)*top1-top2+1;

    return total;
}

__device__ double rastringin (double solution[], int dimensions)
{
    double total = 0;

    for(int i=0;i<dimensions;i++)
    {
        total=total+(pow(solution[i],(double)2)-10*cos(2*M_PI*solution[i])+10);
    }

    return total;
}

__device__ double ackley(double solution[], int dimensions)
{
    int i;
    double aux, aux1, result;

    for (i = 0; i < dimensions; i++)
    {
        aux += solution[i]*solution[i];
    }
    for (i = 0; i < dimensions; i++)
    {
        aux1 += cos(2.0*M_PI*solution[i]);
    }

    result = -20.0*(exp(-0.2*sqrt(1.0/(float)dimensions*aux)))-exp(1.0/(float)dimensions*aux1)+20.0+exp(1.0);

    return result;
}


__device__ double rosenbrock (double solution[], int dimensions)
{
    double total = 0;
    for (int i = 0; i < dimensions-1; i++)
    {
        total=total+100.*pow((solution[i+1] - pow(solution[i],2.)),2) + pow((1. - solution[i]),2);
    }

    return total;
}


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
    for (int i = 0; i < dbats_count; i++) {
        if (bats[i].fitness < current_best_val) {
            current_best_val = bats[i].fitness;
            best_indice = i;
        }
    }
    copy_bat(&bats[best_indice], best);
}


__device__ void log_bat_stdout(struct bat *bat, int dimensions) 
{
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        position_average+=bat->position[i];
    }
    printf("ITERATIONS: %d\n", diterations);
    printf("BATS_COUNT: %d\n", dbats_count);
    printf("DIMENSIONS: %d\n", DIMENSIONS);
    printf("Fitness E: %E\n", bat->fitness);
}


__device__ double  my_rand(curandState *state, double min, double max)
{
    float myrandf = curand_uniform(&state[blockIdx.x]);

    myrandf *= (max - min );
    myrandf += min;

    return myrandf;
}


__device__ int my_rand_int(curandState *state, double min, double max)
{
    float myrandf = curand_uniform(&state[blockIdx.x]);
    myrandf *= (max - min );
    myrandf += min;
    return (int)truncf(myrandf);
}

__device__ void initialize_function(void)
{
    switch(EVALUTAION_FUNCTION) {
        case SPHERE:
            BOUNDRY_MIN = 0.0;
            BOUNDRY_MAX = 100.0;
            objective_function = &sphere;
            break;
        case RASTRINGIN:
            BOUNDRY_MIN = -5.12;
            BOUNDRY_MAX = 5.12;
            objective_function = &rastringin;
            break;
        case GRIEWANK:
            BOUNDRY_MIN = -600.0;
            BOUNDRY_MAX = 600.0;
            objective_function = &griewank;
            break;
        case ACKLEY:
            BOUNDRY_MIN = -32.0;
            BOUNDRY_MAX = 32.0;
            objective_function = &ackley;
            break;
        case ROSENBROOK:
            BOUNDRY_MIN = -30.0;
            BOUNDRY_MAX = 30.0;
            objective_function = &rosenbrock;
            break;
    }
}
__global__ void run_bats(curandState *state, unsigned int seed, struct bat *bats, struct bat *candidates, int iterations, int bats_count)
{
    dbats_count = bats_count;
    diterations = iterations;
    initialize_function();
    int id = threadIdx.x + blockIdx.x * 64;
    curand_init(seed, id, 0, &state[id]);
    curandState localState = state[id];

    __shared__ struct bat *best;
    __shared__ int iteration;
    __shared__ double loudness_average;
    best = (struct bat *) malloc(sizeof(struct bat));

    loudness_average = 1.0;
    bats[threadIdx.x].pulse_rate = 0.0;
    bats[threadIdx.x].frequency = 0.0;
    bats[threadIdx.x].loudness = INITIAL_LOUDNESS;


    for (int j = 0; j < DIMENSIONS; j++) {
        bats[threadIdx.x].velocity[j] = 0;
        bats[threadIdx.x].position[j] = my_rand(&localState, BOUNDRY_MIN, BOUNDRY_MAX);

    }

    bats[threadIdx.x].fitness = objective_function(bats[threadIdx.x].position, DIMENSIONS);

     __syncthreads();

    get_best(bats, best);

    iteration = 0;
    while(iteration < iterations) {
            //frequency
            double beta = my_rand(&localState, BETA_MIN, BETA_MAX);
            beta = FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta;

            bats[threadIdx.x].frequency = beta;
            //velocity
            for (int i = 0; i < DIMENSIONS; ++i) {
                bats[threadIdx.x].velocity[i]+=  (bats[threadIdx.x].position[i] - best->position[i]) * bats[threadIdx.x].frequency;

                if (bats[threadIdx.x].velocity[i] > BOUNDRY_MAX) {
                    bats[threadIdx.x].velocity[i] = BOUNDRY_MAX;
                } else if (bats[threadIdx.x].velocity[i] < BOUNDRY_MIN) {
                    bats[threadIdx.x].velocity[i] = BOUNDRY_MIN;
                }
            }

            copy_bat(&bats[threadIdx.x], &candidates[threadIdx.x]);

            //update position
            for (int i = 0; i < DIMENSIONS; ++i) {
                candidates[threadIdx.x].position[i] += candidates[threadIdx.x].velocity[i];

                if (candidates[threadIdx.x].position[i] > BOUNDRY_MAX) {
                    candidates[threadIdx.x].position[i] = BOUNDRY_MAX;
                } else if (candidates[threadIdx.x].position[i] < BOUNDRY_MIN) {
                    candidates[threadIdx.x].position[i] = BOUNDRY_MIN;
                }

            }

            //local search
            if (my_rand(&localState, 0.0, 1.0) < candidates[threadIdx.x].pulse_rate) {
                for (int i = 0; i < DIMENSIONS; i++ ) {
                    candidates[threadIdx.x].position[i] = best->position[i] +  loudness_average * my_rand(&localState, -1.0, 1.0);
                }
            }

            //position perturbation
            int dimension = my_rand_int(&localState, 0, DIMENSIONS);
            candidates[threadIdx.x].position[dimension] = candidates[threadIdx.x].position[dimension] * my_rand(&localState, 0.0,1.0);


            bats[threadIdx.x].fitness = objective_function(bats[threadIdx.x].position, DIMENSIONS);
            candidates[threadIdx.x].fitness = objective_function(candidates[threadIdx.x].position, DIMENSIONS);

            if (my_rand(&localState, 0.0,1.0) < bats[threadIdx.x].loudness || candidates[threadIdx.x].fitness < bats[threadIdx.x].fitness) {
                copy_bat(&candidates[threadIdx.x], &bats[threadIdx.x]);
                bats[threadIdx.x].pulse_rate = 1 - exp(-LAMBDA*iteration);
                bats[threadIdx.x].loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
            }
            loudness_average=0;
            loudness_average+=bats[threadIdx.x].loudness;
            __syncthreads();
            get_best(bats, best);
            loudness_average/= dbats_count;

            iteration++;
    }

    if (threadIdx.x == 0) {
        log_bat_stdout(best, DIMENSIONS);
    }

    __syncthreads();
}

int main(int argc, char **argv)
{
    clock_t begin = clock();
    char *HELP = "--help";

    if (argc > 1 && strcmp(argv[1], HELP) == 0) {
        printf("The GPU version of the BAT algorithm");
        return 0;
    }

    char* sIterations;
    sIterations = getenv("ITERATIONS");
    if (sIterations != NULL) {
        iterations = atoi(sIterations);
    }

    char* sBatsCount;
    sBatsCount = getenv("BATS_COUNT");
    if (sBatsCount != NULL) {
        bats_count = atoi(sBatsCount);
    }


    struct bat *bats;
    struct bat *candidates;
    int size_of_bats = bats_count * sizeof(struct bat) ;

    CUDA_CALL(cudaMalloc, (void **)&bats, size_of_bats);
    CUDA_CALL(cudaMalloc, (void **)&candidates, size_of_bats);


    curandState *deviceStates;
    CUDA_CALL(cudaMalloc, (void **)&deviceStates, bats_count *sizeof(curandState));

    run_bats<<<1,bats_count>>>(deviceStates, time(NULL), bats, candidates, iterations, bats_count);

    CUDA_CALL(cudaDeviceSynchronize);
    CUDA_CALL(cudaFree, bats);
    CUDA_CALL(cudaFree, candidates);
    CUDA_CALL(cudaFree, deviceStates);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Function %s\n", get_function_name(EVALUTAION_FUNCTION));
    printf("Time took GPU: %f\n", time_spent);

    return 0;
}
