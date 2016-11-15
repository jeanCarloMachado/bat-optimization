#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include "common.h"

#define DIMENSIONS 30
#define MAX_ITERATIONS 1000
//#define MAX_ITERATIONS 1000
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

enum {LOG_OBJECTIVE, LOG_RANDOM, LOG_STDOUT, LOG_SCALAR_ATRIBUTES, LOG_VECTOR_ATRIBUTES};

extern double rosenbrock (double solution[], int dimensions);
extern double sphere (double solution[], int dimensions);
extern double schewefel (double solution[], int dimensions);
extern double ackley (double solution[], int dimensions);
extern double rastringin (double solution[], int dimensions);
extern double griewank (double solution[], int dimensions);

enum {ROSENBROOK, SPHERE, SCHEWEFEL, ACKLEY, RASTRINGIN, GRIEWANK};
double (*objective_function)(double[], int);

int BOUNDRY_MAX;
int BOUNDRY_MIN;

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

void initialize_bats(struct bat *bats)
{

    for (int i = 0; i < BATS_COUNT; i++) {
        bats[i].pulse_rate = 0;
        bats[i].frequency = 0;
        bats[i].fitness = 0;
        bats[i].loudness = INITIAL_LOUDNESS;

        for (int j = 0; j < DIMENSIONS; j++) {
            bats[i].velocity[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
            bats[i].position[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
        }
    }
}

void log_bat(struct bat *bat)
{
    logger(LOG_SCALAR_ATRIBUTES, "F,PR,L: %f %f %f\n", bat->frequency, bat->pulse_rate, bat->loudness);

    for (int i = 0; i < DIMENSIONS; i++) {
        logger(LOG_VECTOR_ATRIBUTES, "%f\t%f\t%f\n", bat->velocity[i], bat->position[i], bat->fitness);
    }
}

struct bat get_best(struct bat *bats, struct bat *best)
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

void update_velocity(struct bat *bat, struct bat *best)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->velocity[i] = bat->velocity[i] + (bat->position[i] - best->position[i]) * bat->frequency;

        if (bat->velocity[i] > BOUNDRY_MAX) {
            bat->velocity[i] = BOUNDRY_MAX;
        }
        if (bat->velocity[i] < BOUNDRY_MIN) {
            bat->velocity[i] = BOUNDRY_MIN;
        }
    }
}

double generate_frequency()
{
    double beta = my_rand(BETA_MIN, BETA_MAX);
    return FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta;
}

void update_position(struct bat *bat)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = bat->position[i] + bat->velocity[i];

        if (bat->position[i] > BOUNDRY_MAX) {
            bat->position[i] = BOUNDRY_MAX;
        }
        if (bat->position[i] < BOUNDRY_MIN) {
            bat->position[i] = BOUNDRY_MIN;
        }
    }
}

void local_search(struct bat *bat, struct bat *best, double loudness_average)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = best->position[i] + loudness_average * my_rand(0.0,1.0);
    }
}

double calc_loudness_average(struct bat *bats)
{
    double total = 0;


    for(int i=0;i<BATS_COUNT;i++) {
        total+= bats[i].loudness;
    }

    return total / BATS_COUNT;
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

void decrease_loudness(struct bat *bat, int iteration)
{
    //linear method
    //bat->loudness = INITIAL_LOUDNESS - ((INITIAL_LOUDNESS/100)*iteration);

    //geometric method
    bat->loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
}

void position_perturbation(struct bat *bat)
{
    int dimension = my_rand(0, DIMENSIONS);
    bat->position[dimension] = bat->position[dimension] * my_rand(0,1);
}

void force_boundry_over_position(struct bat *bat)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        if (bat->position[i] > BOUNDRY_MAX) {
            bat->position[i] = BOUNDRY_MAX;
        }
        if (bat->position[i] < BOUNDRY_MIN) {
            bat->position[i] = BOUNDRY_MIN;
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




int main()
{
    struct bat *bats;
    struct bat *best;
    struct bat *candidate;
    int iteration;
    double best_result,average_result,worst_result;


    allocate_resources();


    const int EVALUTAION_FUNCTION = ROSENBROOK;
    switch(EVALUTAION_FUNCTION) {
        case SPHERE:
            BOUNDRY_MIN = -100;
            BOUNDRY_MAX = 100;
            objective_function = &sphere; 
            break;
        case RASTRINGIN:
            BOUNDRY_MIN = -5.12;
            BOUNDRY_MAX = 5.12;
            objective_function = &rastringin; 
            break;
        case GRIEWANK:
            BOUNDRY_MIN = -600;
            BOUNDRY_MAX = 600;
            objective_function = &griewank; 
            break;
        case ACKLEY:
            BOUNDRY_MIN = -32;
            BOUNDRY_MAX = 32;
            objective_function = &ackley; 
            break;
        case SCHEWEFEL:
            BOUNDRY_MIN = -500;
            BOUNDRY_MAX = 500;
            objective_function = &schewefel; 
            break;
        case ROSENBROOK:
            BOUNDRY_MIN = -30;
            BOUNDRY_MAX = 30;
            objective_function = &rosenbrock; 
            break;
    }


    my_seed();

    bats = (struct bat *) malloc(sizeof(struct bat) * BATS_COUNT);
    best = (struct bat *) malloc(sizeof(struct bat));
    candidate = (struct bat *) malloc(sizeof(struct bat));

    initialize_bats(bats);

    get_best(bats, best); 

    for (iteration = 0; iteration < MAX_ITERATIONS ; ++iteration) {
        for (int j = 0; j < BATS_COUNT; j++) {
            bats[j].frequency = generate_frequency();
            update_velocity(&bats[j], best);
            memcpy(candidate, &bats[j], sizeof(struct bat));

            update_position(candidate);

            if (my_rand(0,1) < candidate->pulse_rate) {
                local_search(candidate, best, calc_loudness_average(bats));
            }

            position_perturbation(candidate);
            force_boundry_over_position(candidate);

            bats[j].fitness = objective_function(bats[j].position, DIMENSIONS);
            candidate->fitness = objective_function(candidate->position, DIMENSIONS);
            if (my_rand(0,1) < bats[j].loudness || candidate->fitness < bats[j].fitness) {
                memcpy(bats[j].position, candidate->position, sizeof candidate->position);
                bats[j].fitness = candidate->fitness;
                bats[j].pulse_rate = 1 - exp(-LAMBDA*iteration);
                decrease_loudness(&bats[j], iteration);
            }
            get_best(bats, best);
            if (LOG_ATRIBUTES_ENABLED) {
                log_bat(&bats[j]);
            }
        }
        get_best(bats, best);

        if (LOG_OBJECTIVE_ENABLED) {
            average_result = fitness_average(bats);
            worst_result = get_worst(bats).fitness;
            logger(
                    LOG_OBJECTIVE,
                    "%f\t%f\t%f\n",
                    best->fitness,
                    average_result,
                    worst_result
                  );

        }

        /* if (fabs(best_result) < 0.0009) { */
        /*     break; */
        /* } */

    }

    logger(
            LOG_STDOUT,
            "Best of All: %f iterations (%d)",
            best->fitness,
            iteration
          );

    free(bats);
    free(best);
    free(candidate);
    deallocate_resources();
    return 0;
}


