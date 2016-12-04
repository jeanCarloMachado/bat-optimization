#include <stdio.h>
#include <unistd.h>
#include <curand.h>
#include <curand_kernel.h>

#define ITERATIONS 10000
#define BATS_COUNT 256
#define DIMENSIONS 100

#define INITIAL_LOUDNESS 1.0

//probability of accepting bad results #define ALFA 0.5
//affects local search
#define LAMBDA 0.1
#define ALFA 0.5

#define BETA_MAX 1.0
#define BETA_MIN 0.0

enum {ROSENBROOK, SPHERE, SCHWEFEL, ACKLEY, RASTRINGIN, GRIEWANK, SHUBER};
const int EVALUTAION_FUNCTION = ROSENBROOK;

__device__ int BOUNDRY_MAX;
__device__ int BOUNDRY_MIN;
__device__ int FREQUENCY_MIN;
__device__ int FREQUENCY_MAX;

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

__device__ double (*objective_function)(double[], int);

__device__ double schwefel(double solution[], int dimensions)
{
    double aux = 0;
    for (int i=0;i<dimensions;i++)
    {
        aux += solution[i]*sin(sqrt(fabs(solution[i]))); 
    }
    return(-1*aux/dimensions);
}


__device__ double shuber (double solution[], int dimensions)
{
    /*
       -   Domain  |x| <= 10.0
       -   Number of local minimum = 400
       -   Global minimum fmin = -24.062499 at the ff. points
       -    (-6.774576, -6.774576), ..., (5.791794, 5.791794)
       */
    double sum = 0.0;
    for (int i = 0; i < dimensions; i++) {
        sum += -sin(2.0*solution[i]+1.0)
            -2.0*sin(3.0*solution[i]+2.0)
            -3.0*sin(4.0*solution[i]+3.0)
            -4.0*sin(5.0*solution[i]+4.0)
            -5.0*sin(6.0*solution[i]+5.0);
    }
    return sum;
}

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

    //result = -20.0*(exp(-0.2*sqrt(1.0/(float)dimensions*aux)))-exp(1.0/(float)dimensions*aux1)+20.0+exp(1);
    result = 666;

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
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        position_average+=bat->position[i];
    }
    /* printf("Frequency: %f\n", bat->frequency); */
    /* printf("Pulse-rate: %f\n", bat->pulse_rate); */
    /* printf("Loudness: %f\n", bat->loudness); */
    printf("ITERATIONS: %d\n", ITERATIONS);
    printf("BATS_COUNT: %d\n", BATS_COUNT);
    printf("DIMENSIONS: %d\n", DIMENSIONS);
    printf("Fitness E: %E\n", bat->fitness);
}



__device__ double  my_rand(curandState_t state, double min, double max)
{
    float myrandf = curand_uniform(&state);
    myrandf *= (max - min );
    myrandf += min;

    return myrandf;
}


__device__ int my_rand_int(curandState_t state, double min, double max)
{
    float myrandf = curand_uniform(&state);
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
        case SHUBER:
            BOUNDRY_MIN = -100.0;
            BOUNDRY_MAX = 100.0;
            objective_function = &shuber;
            break;
        case SCHWEFEL:
            BOUNDRY_MIN = -500.0;
            BOUNDRY_MAX = 500.0;
            objective_function = &schwefel;
            break;
        case ROSENBROOK:
            BOUNDRY_MIN = -30.0;
            BOUNDRY_MAX = 30.0;
            objective_function = &rosenbrock;
            break;
    }
}
__global__ void run_bats(unsigned int seed, struct bat *bats, struct bat *candidates)
{
    initialize_function();
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

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
        bats[threadIdx.x].position[j] = my_rand(state, BOUNDRY_MIN, BOUNDRY_MAX);
    }

    bats[threadIdx.x].fitness = objective_function(bats[threadIdx.x].position, DIMENSIONS);

    __syncthreads();

    get_best(bats, best);

    iteration = 0;
    while(iteration < ITERATIONS) {
            //frequency
            double beta = my_rand(state, BETA_MIN, BETA_MAX);
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
            if (my_rand(state, 0.0, 1.0) < candidates[threadIdx.x].pulse_rate) {
                for (int i = 0; i < DIMENSIONS; i++ ) {
                    candidates[threadIdx.x].position[i] = best->position[i] +  loudness_average * my_rand(state, -1.0, 1.0);
                }
            }


            //position perturbation
            int dimension = my_rand_int(state, 0, DIMENSIONS);
            candidates[threadIdx.x].position[dimension] = candidates[threadIdx.x].position[dimension] * my_rand(state, 0.0,1.0);


            bats[threadIdx.x].fitness = objective_function(bats[threadIdx.x].position, DIMENSIONS);
            candidates[threadIdx.x].fitness = objective_function(candidates[threadIdx.x].position, DIMENSIONS);

            if (my_rand(state, 0.0,1.0) < bats[threadIdx.x].loudness || candidates[threadIdx.x].fitness < bats[threadIdx.x].fitness) {
                /* printf("I: %d\n", iteration); */
                copy_bat(&candidates[threadIdx.x], &bats[threadIdx.x]);
                bats[threadIdx.x].pulse_rate = 1 - exp(-LAMBDA*iteration);
                bats[threadIdx.x].loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
            }
            loudness_average=0;
            loudness_average+=bats[threadIdx.x].loudness;
            get_best(bats, best);
            __syncthreads();
            loudness_average/= BATS_COUNT;

            iteration++;
    }

    if (threadIdx.x == 0) {
        log_bat_stdout(best, DIMENSIONS);
    }
}


int main(void)
{
    struct bat *bats;
    struct bat *candidates;
    cudaMalloc((void **)&bats, sizeof(struct bat) * BATS_COUNT);
    cudaMalloc((void **)&candidates, sizeof(struct bat) * BATS_COUNT);

    run_bats<<<1,BATS_COUNT>>>(time(NULL), bats, candidates);
    cudaFree(bats);

    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
        printf("kernel launch failed with error \"%s\".\n",
               cudaGetErrorString(cudaerr));

    return 0;
}

