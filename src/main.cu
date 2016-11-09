#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
extern "C" {
#include "common.h"
}

//rastringin
/* #define BOUNDRY_MIN -5.12 */
/* #define BOUNDRY_MAX 5.12 */

//griewank
/* #define BOUNDRY_MIN -600 */
/* #define BOUNDRY_MAX 600 */

//sphere
/* #define BOUNDRY_MIN -100 */
/* #define BOUNDRY_MAX 100 */

//ackley
#define BOUNDRY_MIN -32
#define BOUNDRY_MAX 32

#define DIMENSIONS 200
#define MAX_ITERATIONS 1000
#define BATS_COUNT 40
#define FREQUENCY_MIN 0.0
#define FREQUENCY_MAX 1.0
#define INITIAL_LOUDNESS 1.0

#define DUMP_DIR "./dump"
#define BETA_MIN 0.0
#define BETA_MAX 1.0

//probability of accepting bad results
#define ALFA 0.5
//affects local search
#define LAMBDA 0.1

#define LOG_OBJECTIVE_ENABLED 1
#define LOG_ATRIBUTES_ENABLED 1
#define LOG_RANDOM_ENABLED 0

#define LOG_OBJECTIVE 1
#define LOG_RANDOM 2
#define LOG_STDOUT 3
#define LOG_SCALAR_ATRIBUTES 4
#define LOG_VECTOR_ATRIBUTES 5

int RUN_TIME;
FILE *LOG_OBJECTIVE_FILE;
FILE *LOG_SCALAR_ATRIBUTES_FILE;
FILE *LOG_VECTOR_ATRIBUTES_FILE;
FILE *LOG_RANDOM_FILE;

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

void logger(int destination, char *fmt, ...)
{
    char formatted_string[6666];

    va_list argptr;
    va_start(argptr,fmt);
    vsprintf(formatted_string, fmt, argptr);
    va_end(argptr);

    if (destination == LOG_OBJECTIVE) 
        fprintf(LOG_OBJECTIVE_FILE,"%s",formatted_string);
    else if (destination == LOG_SCALAR_ATRIBUTES)
        fprintf(LOG_SCALAR_ATRIBUTES_FILE,"%s",formatted_string);
    else if (destination == LOG_VECTOR_ATRIBUTES)
        fprintf(LOG_VECTOR_ATRIBUTES_FILE,"%s",formatted_string);
    else if (destination == LOG_RANDOM)
        fprintf(LOG_RANDOM_FILE,"%s",formatted_string);
    else if (destination == LOG_STDOUT) 
        printf("%s",formatted_string);
}

void my_rand_collection(int min, int max, double *collection, int size)
{
    for (int i = 0; i < size; ++i) {
        collection[i]  = my_rand(min, max);
    }
}

__global__ void initialize_bat(struct bat *bats, double *velocity, double *position, int size)
{
    bats[blockIdx.x].pulse_rate = 0;
    bats[blockIdx.x].frequency = 0;
    bats[blockIdx.x].fitness = 0;
    bats[blockIdx.x].loudness = INITIAL_LOUDNESS;

    memcpy(bats[blockIdx.x].velocity, &velocity[blockIdx.x * DIMENSIONS], DIMENSIONS * sizeof(double));
    memcpy(bats[blockIdx.x].position, &position[blockIdx.x * DIMENSIONS], DIMENSIONS * sizeof(double));
}

void initialize_bats(struct bat *bats)
{
    double *velocity, *position;
    double *d_velocity, *d_position;

    int vector_positions,size;
    vector_positions = DIMENSIONS * BATS_COUNT;
    size = sizeof(double) * vector_positions;

    velocity  = (double *)malloc(size);
    position  = (double *)malloc(size);

    cudaMalloc((void **)&d_velocity, size);
    cudaMalloc((void **)&d_position, size);

    my_rand_collection(BOUNDRY_MIN, BOUNDRY_MAX, velocity, vector_positions);
    my_rand_collection(BOUNDRY_MIN, BOUNDRY_MAX, position, vector_positions);

    cudaMemcpy(d_velocity, velocity, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_position, position, size, cudaMemcpyHostToDevice);

    initialize_bat<<<BATS_COUNT,1>>>(bats, d_velocity, d_position, vector_positions);

    free(velocity);
    free(position);
    cudaFree(d_velocity);
    cudaFree(d_position);
}
void log_bat(struct bat *bat)
{
    struct bat *l_bat; 
    l_bat = (struct bat *) malloc(sizeof(struct bat));
    cudaMemcpy(l_bat, bat, sizeof(struct bat), cudaMemcpyDeviceToHost);


    logger(LOG_SCALAR_ATRIBUTES, "F,PR,L: %f %f %f\n", l_bat->frequency, l_bat->pulse_rate, l_bat->loudness);

    for (int i = 0; i < DIMENSIONS; i++) {
        logger(LOG_VECTOR_ATRIBUTES, "%f\t%f\t%f\n", l_bat->velocity[i], l_bat->position[i], l_bat->fitness);
    }

    free(l_bat);
}

__global__ void get_best(struct bat *bats, struct bat *best)
{
    double current_best_val; 
    double current_val;

    current_val = current_best_val = bats[0].fitness;
    memcpy(best, &bats[0], sizeof(struct bat));
    for (int i = 0; i < BATS_COUNT; i++) {
        current_val = bats[i].fitness;
        if (current_val < current_best_val) {
            current_best_val = current_val;
            memcpy(best, &bats[i], sizeof(struct bat));
        }
    }
}

__global__ void update_velocity(struct bat *bat, struct bat *best, double random)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->velocity[i] = bat->velocity[i] + (bat->position[i] - best->position[i]) * bat->frequency;
        if (bat->velocity[i] > BOUNDRY_MAX || bat->velocity[i] < BOUNDRY_MIN) {
            bat->velocity[i] = random;
        }
    }
}

__global__ void generate_frequency(struct bat *bat, double beta)
{
    bat->frequency = FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta;
}

__global__ void update_position(struct bat *bat, double *random_collection)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = bat->position[i] + bat->velocity[i];

        if (bat->position[i] > BOUNDRY_MAX || bat->position[i] < BOUNDRY_MIN) {
            bat->position[i] = random_collection[i];
        }
    }
}

__device__ float calc_loudness_average(struct bat *bats)
{
    double total = 0;


    for(int i=0;i<BATS_COUNT;i++) {
        total+= bats[i].loudness;
    }

    return total / BATS_COUNT;
}


__global__ void local_search(struct bat *bat, struct bat *best, struct bat *candidate, struct bat *bats, double rand, double *random_collection)
{
    if (rand >= candidate->pulse_rate) {
        return;
    }

    double loudness_average = calc_loudness_average(bats);
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = best->position[i] + loudness_average * random_collection[i];
    }
}

double fitness_average(struct bat bats[])
{
    double result = 0;

    for (int i = 0; i < BATS_COUNT; i++) {
        result+= bats[i].fitness;
    }

    return result / BATS_COUNT;
}


void allocate_resources()
{
    RUN_TIME = time(NULL);
    char fileName[100];

    if (LOG_OBJECTIVE_ENABLED) {
        sprintf(fileName, "%s/%i-objective", DUMP_DIR, RUN_TIME);
        LOG_OBJECTIVE_FILE = fopen(fileName,"w");
        if (LOG_OBJECTIVE_FILE == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Objective log: %s\n", fileName);
    }

    if (LOG_ATRIBUTES_ENABLED) {
        sprintf(fileName, "%s/%i-scalar_attr", DUMP_DIR, RUN_TIME);
        LOG_SCALAR_ATRIBUTES_FILE = fopen(fileName,"w");
        if (LOG_SCALAR_ATRIBUTES_FILE == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Scalar atributes log: %s\n", fileName);


        sprintf(fileName, "%s/%i-vector_attr", DUMP_DIR, RUN_TIME);
        LOG_VECTOR_ATRIBUTES_FILE = fopen(fileName,"w");
        if (LOG_VECTOR_ATRIBUTES_FILE == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Vector Atributes log: %s\n", fileName);

    }

    if (LOG_RANDOM_ENABLED) {
        sprintf(fileName, "%s/%i-random", DUMP_DIR, RUN_TIME);
        LOG_RANDOM_FILE = fopen(fileName,"w");
        if (LOG_RANDOM_FILE == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Random log: %s\n", fileName);
    }
}

void deallocate_resources()
{
    if (LOG_OBJECTIVE_ENABLED) {
        fclose(LOG_OBJECTIVE_FILE);
    }
    if (LOG_ATRIBUTES_ENABLED) {
        fclose(LOG_SCALAR_ATRIBUTES_FILE);
        fclose(LOG_VECTOR_ATRIBUTES_FILE);
    }
    if (LOG_RANDOM_ENABLED) {
        fclose(LOG_RANDOM_FILE);
    }
}

__device__ void decrease_loudness(struct bat *bat, int iteration)
{
    //linear method
    //bat->loudness = INITIAL_LOUDNESS - ((INITIAL_LOUDNESS/100)*iteration);

    //geometric method
    bat->loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
}

__global__  void position_perturbation(struct bat *bat, int dimension, double random)
{
    bat->position[dimension] = bat->position[dimension] * random;
}

__global__ void force_boundry_over_position(struct bat *bat, double *random_collection)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        if (bat->position[i] > BOUNDRY_MAX || bat->position[i] < BOUNDRY_MIN) {
            bat->position[i] = random_collection[i];
        }
    }
}

struct bat get_worst(struct bat bats[])
{
    double current_worst_val;
    double current_val;

    current_val = current_worst_val = bats[0].fitness;
    struct bat current_worst_bat = bats[0];
    for (int i = 0; i < BATS_COUNT; i++) {
        current_val = bats[i].fitness;
        if (current_worst_val <  current_val) {
            current_worst_val = current_val;
            current_worst_bat = bats[i];
        }
    }

    return current_worst_bat;
}


void objective_function (struct bat *bat)
{

    struct bat *local_bat;

    local_bat = (struct bat *) malloc(sizeof(struct bat));

    cudaMemcpy(local_bat, bat, sizeof(struct bat), cudaMemcpyDeviceToHost);

    //bat->fitness = rastringin(bat->position, DIMENSIONS);
    local_bat->fitness = griewank(local_bat->position, DIMENSIONS);
    /* bat->fitness = sphere(bat->position, DIMENSIONS); */

    //bat->fitness = ackley(bat->position, DIMENSIONS);
    //usleep(0);

    /* printf("Fitness: %f\n", local_bat->fitness); */

    cudaMemcpy(bat, local_bat, sizeof(struct bat), cudaMemcpyHostToDevice);
    free(local_bat);
}


__global__ void greedy_update_on_bat(struct bat *bats, struct bat *candidate, double random, int iteration)
{
    /* printf("Greedy Fitness: %f\n", bats[iteration].fitness); */
    /* printf("Candidate Fitness: %f\n", candidate->fitness); */
    if (random < bats[iteration].loudness || candidate->fitness < bats[iteration].fitness) {
        memcpy(bats[iteration].position, candidate->position, sizeof candidate->position);
        bats[iteration].fitness = candidate->fitness;
        bats[iteration].pulse_rate = 1 - exp(-LAMBDA*iteration);
        decrease_loudness(&bats[iteration], iteration);
    }
}

int main()
{
    allocate_resources();
    struct bat *bats;
    struct bat *best;
    struct bat *best_local;
    struct bat *candidate;

    int iteration;
    double best_result,average_result,worst_result;
    double *random_collection, *random_collection_device;

    my_seed();

    best_local = (struct bat *) malloc(sizeof(struct bat));
    cudaMalloc((void **)&bats, sizeof(struct bat) * BATS_COUNT);
    cudaMalloc((void **)&best, sizeof(struct bat));
    cudaMalloc((void **)&candidate, sizeof(struct bat));

    random_collection = (double *)malloc(sizeof(double) * DIMENSIONS);
    cudaMalloc((void **)&random_collection_device, sizeof(double) * DIMENSIONS);
    initialize_bats(bats);

    /* log_bat(&bats[0]); */
    get_best<<<1,1>>>(bats, best); 

    for (iteration = 0; iteration < MAX_ITERATIONS ; ++iteration) {
        for (int j = 0; j < BATS_COUNT; j++) {
            generate_frequency<<<1,1>>>(&bats[j], my_rand(BETA_MIN, BETA_MAX));
            update_velocity<<<1,1>>>(&bats[j], best, my_rand(BOUNDRY_MIN, BOUNDRY_MAX));
            cudaMemcpy(candidate, &bats[j], sizeof(struct bat), cudaMemcpyHostToHost);



            //////////////


            my_rand_collection(BOUNDRY_MIN, BOUNDRY_MAX, random_collection, DIMENSIONS);
            cudaMemcpy(random_collection_device, random_collection, sizeof(double) * DIMENSIONS, cudaMemcpyHostToDevice);

            update_position<<<1,1>>>(candidate, random_collection_device);


            //////////////



            my_rand_collection(0.0, 1.0, random_collection, DIMENSIONS);

            cudaMemcpy(random_collection_device, random_collection, sizeof(double) * DIMENSIONS, cudaMemcpyHostToDevice);


            local_search<<<1,1>>>(candidate, best, candidate, bats, my_rand(0,1), random_collection_device);

            position_perturbation<<<1,1>>>(candidate, my_rand(0, DIMENSIONS), my_rand(0,1));


            //////////////


            my_rand_collection(BOUNDRY_MIN, BOUNDRY_MAX, random_collection, DIMENSIONS);
            cudaMemcpy(random_collection_device, random_collection, sizeof(double) * DIMENSIONS, cudaMemcpyHostToDevice);

            force_boundry_over_position<<<1,1>>>(candidate, random_collection_device);


            objective_function(&bats[j]);
            objective_function(candidate);

            greedy_update_on_bat<<<1,1>>>(bats, candidate, my_rand(0,1), iteration);
            /* if (LOG_ATRIBUTES_ENABLED) { */
            /*     log_bat(&bats[j]); */
            /* } */
    /*     best_result = best->fitness; */

            cudaMemcpy(best_local, best, sizeof(struct bat), cudaMemcpyDeviceToHost);
            if (LOG_OBJECTIVE_ENABLED) {
                /* average_result = fitness_average(bats); */
                /* worst_result = get_worst(bats).fitness; */
                logger(LOG_OBJECTIVE, "%f\n", best_local->fitness);

            }
        }

        /* if (fabs(best_result) < 0.0009) { */
        /*     break; */
        /* } */

    }

    logger(
            LOG_STDOUT,
            "Best of All: %f iterations (%d)",
            best_local->fitness,
            iteration
          );

    cudaFree(bats);
    cudaFree(best);
    cudaFree(candidate);
    cudaFree(random_collection_device);
    free(random_collection);
    free(best_local);
    deallocate_resources();
    return 0;
}


