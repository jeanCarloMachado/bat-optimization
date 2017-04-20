#include <unistd.h>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "internal.h"

#define LAMBDA 0.1
#define ALFA 0.5
#define BETA_MAX 1.0
#define BETA_MIN 0.0
#define INITIAL_LOUDNESS 1.0


extern int dimensions;
extern int bats_count;
extern int evaluation_function;
extern int iterations;

__device__ int devaluation_function;
__device__ int dbats_count;
__device__ int diterations;
enum { N = 624 };        // length of state vector
enum { M = 397 };        // period parameter
__device__ unsigned long state[N];  // internal state
__device__ unsigned long *pNext;    // next value to get from state
__device__ int left;          	 // number of values left before reload needed
__device__ unsigned long MT_randInt(unsigned long n);
__device__ unsigned long randInt();
__device__ void reload();
__device__ unsigned long twist(unsigned long m, unsigned long s0, unsigned long s1);
__device__ unsigned long hiBit(unsigned long u);
__device__ unsigned long loBit(unsigned long u);
__device__ unsigned long loBits(unsigned long u);
__device__ unsigned long mixBits(unsigned long u, unsigned long v );
__device__ void MT_seed(time_t time, clock_t clock);
__device__ unsigned long MT_hash(time_t t, clock_t c);
__device__ void MT_seedfinal(unsigned long oneSeed);
__device__ void MT_initialize(unsigned long seed);
__device__ float MT_randfloat();
__device__ double MT_randExc(const double  *n );

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

typedef struct Bat {
    double pulse_rate;
    double loudness;
    double fitness;
    double frequency;
    double position[1000];
    double velocity[1000];
} Bat;

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

__device__ void copy_bat(struct Bat *from, struct Bat *to)
{
    memcpy(to, from, sizeof(struct Bat));
}

__device__ void  get_best(struct Bat *bats, struct Bat *best)
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


__device__ void log_bat_stdout(struct Bat *bat, int dimensions) 
{
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        position_average+=bat->position[i];
    }
    printf("ITERATIONS: %d\n", diterations);
    printf("BATS_COUNT: %d\n", dbats_count);
    printf("DIMENSIONS: %d\n", dimensions);
    printf("POPULATION: %d\n", dbats_count);
    printf("Fitness E: %E\n", bat->fitness);
}


__device__ double  my_rand(double inferior, double superior)

{
  double result = (double)inferior + ((superior - inferior)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));

    return result;
}


__device__ void initialize_function(void)
{
    switch(devaluation_function) {
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

__global__ void local_bat_run(time_t time, clock_t clock, struct Bat *bats, struct Bat *candidates, int iterations, int bats_count, int evaluation_function, int dimensions)
{
    MT_seed(time, clock);
    dbats_count = bats_count;
    devaluation_function = evaluation_function;
    diterations = iterations;
    initialize_function();

    __shared__ struct Bat *best;
    __shared__ int iteration;
    __shared__ double loudness_average;
    best = (struct Bat *) malloc(sizeof(struct Bat));

    loudness_average = 1.0;
    bats[threadIdx.x].pulse_rate = 0.0;
    bats[threadIdx.x].frequency = 0.0;
    bats[threadIdx.x].loudness = INITIAL_LOUDNESS;


    for (int j = 0; j < dimensions; j++) {
        bats[threadIdx.x].velocity[j] = 0;
        bats[threadIdx.x].position[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);

    }

    bats[threadIdx.x].fitness = objective_function(bats[threadIdx.x].position, dimensions);

     __syncthreads();

    get_best(bats, best);

    iteration = 0;
    while(iteration < iterations) {
            //frequency
            double beta = my_rand(BETA_MIN, BETA_MAX);

            bats[threadIdx.x].frequency = FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta;

            //velocity
            for (int i = 0; i < dimensions; ++i) {
                bats[threadIdx.x].velocity[i]+=  (bats[threadIdx.x].position[i] - best->position[i]) * bats[threadIdx.x].frequency;

                if (bats[threadIdx.x].velocity[i] > BOUNDRY_MAX) {
                    bats[threadIdx.x].velocity[i] = BOUNDRY_MAX;
                } else if (bats[threadIdx.x].velocity[i] < BOUNDRY_MIN) {
                    bats[threadIdx.x].velocity[i] = BOUNDRY_MIN;
                }
            }

            copy_bat(&bats[threadIdx.x], &candidates[threadIdx.x]);

            //update position
            for (int i = 0; i < dimensions; ++i) {
                candidates[threadIdx.x].position[i] += candidates[threadIdx.x].velocity[i];

                if (candidates[threadIdx.x].position[i] > BOUNDRY_MAX) {
                    candidates[threadIdx.x].position[i] = BOUNDRY_MAX;
                } else if (candidates[threadIdx.x].position[i] < BOUNDRY_MIN) {
                    candidates[threadIdx.x].position[i] = BOUNDRY_MIN;
                }

            }

            //local search
            if (my_rand(0.0, 1.0) < candidates[threadIdx.x].pulse_rate) {
                for (int i = 0; i < dimensions; i++ ) {
                    candidates[threadIdx.x].position[i] = best->position[i] +  loudness_average * my_rand(-1.0, 1.0);
                }
            }

            //position perturbation
            int dimension = my_rand(0, dimensions);
            candidates[threadIdx.x].position[dimension] = candidates[threadIdx.x].position[dimension] * my_rand(0.0,1.0);


            bats[threadIdx.x].fitness = objective_function(bats[threadIdx.x].position, dimensions);
            candidates[threadIdx.x].fitness = objective_function(candidates[threadIdx.x].position, dimensions);

            if (my_rand(0.0,1.0) < bats[threadIdx.x].loudness && candidates[threadIdx.x].fitness < bats[threadIdx.x].fitness) {
                copy_bat(&candidates[threadIdx.x], &bats[threadIdx.x]);
                bats[threadIdx.x].pulse_rate = 1 - exp(-LAMBDA*iteration);
                bats[threadIdx.x].loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
            }
            loudness_average=0;
            loudness_average+=bats[threadIdx.x].loudness;
            get_best(bats, best);
            __syncthreads();
            loudness_average/= dbats_count;

            iteration++;
    }

    if (threadIdx.x == 0) {
        log_bat_stdout(best, dimensions);
    }

    __syncthreads();
}

void bat_run(void)
{
    struct Bat *bats;
    struct Bat *candidates;
    int size_of_bats = bats_count * sizeof(struct Bat) ;

    CUDA_CALL(cudaMalloc, (void **)&bats, size_of_bats);
    CUDA_CALL(cudaMalloc, (void **)&candidates, size_of_bats);

    curandState *deviceStates;
    CUDA_CALL(cudaMalloc, (void **)&deviceStates, bats_count *sizeof(curandState));

    local_bat_run<<<1,bats_count>>>(time(NULL), clock(), bats, candidates, iterations, bats_count, evaluation_function, dimensions);

    CUDA_CALL(cudaDeviceSynchronize);
    CUDA_CALL(cudaFree, bats);
    CUDA_CALL(cudaFree, candidates);
    CUDA_CALL(cudaFree, deviceStates);
}



__device__ unsigned long MT_randInt(unsigned long n)
{
	unsigned long used = n;
	used |= used >> 1;
	used |= used >> 2;
	used |= used >> 4;
	used |= used >> 8;
	used |= used >> 16;
	unsigned long i;
	do{
		i = randInt() & used;
	}while( i > n );
	return i;
}
__device__ unsigned long randInt()
{
	register unsigned long s1;

	if( left == 0 ) reload();
	--left;

	s1 = *pNext++;
	s1 ^= (s1 >> 11);
	s1 ^= (s1 <<  7) & 0x9d2c5680UL;
	s1 ^= (s1 << 15) & 0xefc60000UL;

	return ( s1 ^ (s1 >> 18) );
}

__device__ void reload()
{
	register unsigned long *p = state;
	register int i;
	for( i = N - M; i--; ++p )
		*p = twist( p[M], p[0], p[1] );
	for( i = M; --i; ++p )
		*p = twist( p[M-N], p[0], p[1] );
	*p = twist( p[M-N], p[0], state[0] );

	left = N, pNext = state;
}
__device__ unsigned long twist(unsigned long m, unsigned long s0, unsigned long s1 )
{
	return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfUL);
}

__device__ void MT_seed(time_t time, clock_t clock)
{
	MT_seedfinal( MT_hash( time, clock ) );
}

__device__ unsigned long MT_hash(time_t t, clock_t c)
{
	size_t i, j;
	static unsigned long  differ = 0;

	unsigned long  h1 = 0;
	unsigned char *p = (unsigned char *) &t;
	for(i = 0; i < sizeof(t); ++i)

	{
		h1 *= UCHAR_MAX + 2U;
		h1 += p[i];
	}
	unsigned long  h2 = 0;
	p = (unsigned char *) &c;
	for(j = 0; j < sizeof(c); ++j)
	{
		h2 *= UCHAR_MAX + 2U;
		h2 += p[j];
	}
	return ( h1 + differ++ ) ^ h2;
}

__device__ void MT_seedfinal(unsigned long oneSeed)
{
	MT_initialize(oneSeed);
	reload();
}

__device__ void MT_initialize(unsigned long seed)
{
	register unsigned long *s = state;
	register unsigned long *r = state;
	register int i = 1;
	*s++ = seed & 0xffffffffUL;
	for( ; i < N; ++i )
	{
		*s++ = ( 1812433253UL * ( *r ^ (*r >> 30) ) + i ) & 0xffffffffUL;
		r++;
	}
}

__device__ float MT_randfloat()
{
	return (float)(randInt()) * (1.0/4294967295.0);
}

__device__ double MT_rand()
    { return (double) (randInt()) * (1.0/4294967296.0);
    }

__device__ double MT_randExc(const double  *n )
    { return MT_rand() * *n;
    }

__device__ unsigned long hiBit(unsigned long u) { return u & 0x80000000UL; }
__device__ unsigned long loBit(unsigned long u) { return u & 0x00000001UL; }
__device__ unsigned long loBits(unsigned long u){ return u & 0x7fffffffUL; }
__device__ unsigned long mixBits(unsigned long u, unsigned long v ) { return hiBit(u) | loBits(v); }
