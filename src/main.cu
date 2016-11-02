#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
extern "C" {
#include "mersenne.h"
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
#define LOG_ATRIBUTES_ENABLED 0
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

struct bat get_worst(struct bat bats[]);
void initialize_bats(struct bat bats[]);
double my_rand(int, int);
void my_rand_collection(int min, int max, double *collection, int size);
void my_seed(void);
void log_bat(struct bat *bat);
struct bat get_best(struct bat bats[]);
void update_velocity(struct bat *bat, struct bat *best);
double generate_frequency();
void update_position(struct bat *bat);
void local_search(struct bat *bat, struct bat *best, double loudness_average);
double calc_loudness_average(struct bat bats[]);
double fitness_average(struct bat bats[]);
void logger(int destination, char *fmt, ...);
void allocate_resources(void);
void deallocate_resources();
void decrease_loudness(struct bat*, int);
void position_perturbation(struct bat *bat);
void force_boundry_over_position(struct bat *bat);

void objective_function (struct bat *bat);

int main()
{
    allocate_resources();
    struct bat bats[BATS_COUNT];
    struct bat best;
    struct bat candidate;
    int iteration;
    double best_result,average_result,worst_result;

    my_seed();

    initialize_bats(bats);
    best = get_best(bats);

    for (iteration = 0; iteration < MAX_ITERATIONS ; ++iteration) {
        for (int j = 0; j < BATS_COUNT; j++) {
            bats[j].frequency = generate_frequency();
            update_velocity(&bats[j], &best);
            candidate = bats[j];

            update_position(&candidate);

            if (my_rand(0,1) < candidate.pulse_rate) {
                local_search(&candidate, &best, calc_loudness_average(bats));
            }

            position_perturbation(&candidate);
            force_boundry_over_position(&candidate);


            objective_function(&bats[j]);
            objective_function(&candidate);

            if (my_rand(0,1) < bats[j].loudness || candidate.fitness < bats[j].fitness) {
                memcpy(bats[j].position, candidate.position, sizeof candidate.position);
                bats[j].fitness = candidate.fitness;
                bats[j].pulse_rate = 1 - exp(-LAMBDA*iteration);
                decrease_loudness(&bats[j], iteration);
            }
            best = get_best(bats);
            if (LOG_ATRIBUTES_ENABLED) {
                log_bat(&bats[j]);
            }
        }

        best_result = get_best(bats).fitness;

        if (LOG_OBJECTIVE_ENABLED) {
            average_result = fitness_average(bats);
            worst_result = get_worst(bats).fitness;
            logger(
                    LOG_OBJECTIVE,
                    "%f\t%f\t%f\n",
                    best_result,
                    average_result,
                    worst_result
                  );

        }

        if (fabs(best_result) < 0.0009) {
            break;
        }

    }

    logger(
            LOG_STDOUT,
            "Best of All: %f iterations (%d)",
            best.fitness,
            iteration
          );

    deallocate_resources();
    return 0;
}

void allocate_resources()
{
    RUN_TIME = time(NULL);
    char fileName[100];

    if (LOG_OBJECTIVE_ENABLED) {
        sprintf(fileName, "%s/%i-objective", DUMP_DIR, RUN_TIME);
        LOG_OBJECTIVE_FILE = fopen(fileName,"w+");
        if (LOG_OBJECTIVE_FILE == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Objective log: %s\n", fileName);
    }

    if (LOG_ATRIBUTES_ENABLED) {
        sprintf(fileName, "%s/%i-scalar_attr", DUMP_DIR, RUN_TIME);
        LOG_SCALAR_ATRIBUTES_FILE = fopen(fileName,"w+");
        if (LOG_SCALAR_ATRIBUTES_FILE == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Scalar atributes log: %s\n", fileName);


        sprintf(fileName, "%s/%i-vector_attr", DUMP_DIR, RUN_TIME);
        LOG_VECTOR_ATRIBUTES_FILE = fopen(fileName,"w+");
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

__global__ void initialize_bat(struct bat *bats, double *velocity, double *position, int size)
{
    bats[blockIdx.x].pulse_rate = 0;
    bats[blockIdx.x].frequency = 0;
    bats[blockIdx.x].fitness = 0;
    bats[blockIdx.x].loudness = INITIAL_LOUDNESS;

     memcpy(bats[blockIdx.x].velocity, &velocity[blockIdx.x * DIMENSIONS], DIMENSIONS * sizeof(double));
     memcpy(bats[blockIdx.x].position, &position[blockIdx.x * DIMENSIONS], DIMENSIONS * sizeof(double));
}


void initialize_bats(struct bat bats[])
{
    bat *cuda_bats;
    int bats_size = sizeof(bat) * BATS_COUNT;
    cudaMalloc((void **)&cuda_bats, bats_size);
    cudaMemcpy(cuda_bats, bats, bats_size, cudaMemcpyHostToDevice);

    double *velocity, *position;
    double *d_velocity, *d_position;
    int total,size;
    total = DIMENSIONS * BATS_COUNT;
    size = sizeof(double) * total;

    velocity  = (double *)malloc(size);
    position  = (double *)malloc(size);

    cudaMalloc((void **)&d_velocity, size);
    cudaMalloc((void **)&d_position, size);

    my_rand_collection(BOUNDRY_MIN, BOUNDRY_MAX, velocity, total);
    my_rand_collection(BOUNDRY_MIN, BOUNDRY_MAX, position, total);

    cudaMemcpy(d_velocity, velocity, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_position, position, size, cudaMemcpyHostToDevice);


    initialize_bat<<<BATS_COUNT,1>>>(cuda_bats, d_velocity, d_position, total);

    cudaMemcpy(bats, cuda_bats, bats_size, cudaMemcpyDeviceToHost);

    free(velocity);
    free(position);
    cudaFree(cuda_bats);
    cudaFree(d_velocity);
    cudaFree(d_position);
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

void position_perturbation(struct bat *bat)
{
    int dimension = my_rand(0, DIMENSIONS);
    bat->position[dimension] = bat->position[dimension] * my_rand(0,1);
}


void local_search(struct bat *bat, struct bat *best, double loudness_average)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = best->position[i] + loudness_average * my_rand(0.0,1.0);
    }
}

double calc_loudness_average(struct bat bats[])
{
    double total = 0;

    for(int i=0;i<BATS_COUNT;i++) {
        total+= bats[i].loudness;
    }

    return total / BATS_COUNT;
}


void update_velocity(struct bat *bat, struct bat *best)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->velocity[i] = bat->velocity[i] + (bat->position[i] - best->position[i]) * bat->frequency;
        if (bat->velocity[i] > BOUNDRY_MAX || bat->velocity[i] < BOUNDRY_MIN) {
            bat->velocity[i] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
        }
    }
}

void decrease_loudness(struct bat *bat, int iteration)
{
    //linear method
    //bat->loudness = INITIAL_LOUDNESS - ((INITIAL_LOUDNESS/100)*iteration);

    //geometric method
    bat->loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
}

void update_position(struct bat *bat)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = bat->position[i] + bat->velocity[i];

        if (bat->position[i] > BOUNDRY_MAX || bat->position[i] < BOUNDRY_MIN) {
            bat->position[i] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
        }
    }
}


void force_boundry_over_position(struct bat *bat)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        if (bat->position[i] > BOUNDRY_MAX || bat->position[i] < BOUNDRY_MIN) {
            bat->position[i] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
        }
    }
}



double generate_frequency()
{
    double beta = my_rand(BETA_MIN, BETA_MAX);
    return FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta;
}

void log_bat(struct bat *bat)
{
    logger(LOG_SCALAR_ATRIBUTES, "F,PR,L: %f %f %f\n", bat->frequency, bat->pulse_rate, bat->loudness);

    for (int i = 0; i < DIMENSIONS; i++) {
        logger(LOG_VECTOR_ATRIBUTES, "%f\t%f\t%f\n", bat->velocity[i], bat->position[i], bat->fitness);
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

struct bat get_best(struct bat bats[])
{
    double current_best_val; 
    double current_val;

    current_val = current_best_val = bats[0].fitness;
    struct bat current_best_bat = bats[0];
    for (int i = 0; i < BATS_COUNT; i++) {
        current_val = bats[i].fitness;
        if (current_val < current_best_val) {
            current_best_val = current_val;
            current_best_bat = bats[i];
        }
    }

    return current_best_bat;
}

void my_seed(void)
{
    MT_seed();
}

double my_rand(int min, int max)
{

    double result = (double)min + ((max - min)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));

    if (LOG_RANDOM_ENABLED) {
        logger(LOG_RANDOM, "%i-%i: %f\n", min, max, result); 
    }

    return result;
}

void my_rand_collection(int min, int max, double *collection, int size)
{
    for (int i = 0; i < size; ++i) {
        collection[i]  = my_rand(min, max);
    }
}

void objective_function (struct bat *bat)
{
    //bat->fitness = rastringin(bat->position, DIMENSIONS);
    bat->fitness = griewank(bat->position, DIMENSIONS);
    /* bat->fitness = sphere(bat->position, DIMENSIONS); */

    //bat->fitness = ackley(bat->position, DIMENSIONS);
    //usleep(0);
}


