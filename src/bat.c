#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include "bat.h"

#define MAX_ITERATIONS 100
#define BATS_COUNT 40
#define INITIAL_LOUDNESS 1.0
#define DIMENSIONS 100

//probability of accepting bad results
#define ALFA 0.9
//affects local search
#define LAMBDA 0.9

#define BETA_MAX 1.0
#define BETA_MIN -1.0

const int LOG_OBJECTIVE_ENABLED=1;
const int LOG_ATRIBUTES_ENABLED=1;
const int LOG_RANDOM_ENABLED=0;

void position_perturbation(struct bat *bat);
void decrease_loudness(struct bat *bat, int iteration);
double fitness_average(struct bat bats[]);
double calc_loudness_average(struct bat *bats);
void local_search(struct bat *bat, struct bat *best, double loudness_average); 
void update_position(struct bat *bat);
double generate_frequency();
void update_velocity(struct bat *bat, struct bat *best);
void force_boundry_on_vector(double vector[]);
void force_boundry_on_value(double* value);
void log_bat(struct bat *bat);
void initialize_bats(struct bat *bats, struct bat *best, struct bat *candidate);
void deallocate_bats(struct bat *bats, struct bat *best, struct bat *candidate);
void copy_bat(struct bat *from, struct bat *to);

double (*objective_function)(double[], int);

int BOUNDRY_MAX;
int BOUNDRY_MIN;
int FREQUENCY_MIN;
int FREQUENCY_MAX;

int BOUNDRY_SCAPE_COUNT = 0;
int BOUNDRY_COUNT = 0;

int run_bats(void)
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
            BOUNDRY_MIN = -10.00;
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

    FREQUENCY_MIN=BOUNDRY_MIN;
    FREQUENCY_MAX=BOUNDRY_MAX;

    my_seed();

    bats = (struct bat *) malloc(sizeof(struct bat) * BATS_COUNT);
    best = (struct bat *) malloc(sizeof(struct bat));
    candidate = (struct bat *) malloc(sizeof(struct bat));
    initialize_bats(bats, best, candidate);

    get_best(bats, best); 

    for (iteration = 0; iteration < MAX_ITERATIONS ; ++iteration) {
        for (int j = 0; j < BATS_COUNT; ++j){
            bats[j].frequency = generate_frequency();
            update_velocity(&bats[j], best);
            copy_bat(&bats[j], candidate);

            update_position(candidate);

            if (my_rand(0.0,1.0) < candidate->pulse_rate) {
                local_search(candidate, best, calc_loudness_average(bats));
            }

            position_perturbation(candidate);

            bats[j].fitness = fabs((double)objective_function(bats[j].position, DIMENSIONS));
            candidate->fitness = fabs((double)objective_function(candidate->position, DIMENSIONS));
            if (my_rand(0.0,1.0) < bats[j].loudness || candidate->fitness < bats[j].fitness) {
                memcpy(bats[j].position, candidate->position, (sizeof(double) * DIMENSIONS));
                bats[j].fitness = candidate->fitness;
                bats[j].pulse_rate = 1 - exp(-LAMBDA*iteration);

                decrease_loudness(&bats[j], iteration);
            }
            get_best(bats, best);
            if (LOG_ATRIBUTES_ENABLED) {
                log_bat(&bats[j]);
            }
        }
        /* get_best(bats, best); */

        if (LOG_OBJECTIVE_ENABLED) {
            average_result = fitness_average(bats);
            worst_result = get_worst(bats)->fitness;
            logger(
                    LOG_OBJECTIVE,
                    "%E\t%E\t%E\n",
                    best->fitness,
                    average_result,
                    worst_result
                  );

        }
    }

    log_bat_stdout(best, DIMENSIONS);
    int percentage = (BOUNDRY_SCAPE_COUNT * 100 / BOUNDRY_COUNT);
    printf("Boundry total: %d,escaped: %d, escaped percentage: %d",BOUNDRY_COUNT, BOUNDRY_SCAPE_COUNT, percentage);

    deallocate_bats(bats, best, candidate);
    deallocate_resources();
    return 0;
}


void initialize_bats(struct bat *bats, struct bat *best, struct bat *candidate)
{

    for (int i = 0; i < BATS_COUNT; i++) {
        bats[i].pulse_rate = 0;
        bats[i].frequency = 0;
        bats[i].fitness = 0;
        bats[i].loudness = INITIAL_LOUDNESS;

        bats[i].velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
        bats[i].position = (double *) malloc(sizeof(double) * DIMENSIONS);
        for (int j = 0; j < DIMENSIONS; j++) {
            bats[i].velocity[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
            bats[i].position[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
        }
    }

    best->velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
    best->position = (double *) malloc(sizeof(double) * DIMENSIONS);

    candidate->velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
    candidate->position = (double *) malloc(sizeof(double) * DIMENSIONS);


    for (int j = 0; j < BATS_COUNT; j++) {
        bats[j].fitness = fabs(objective_function(bats[j].position, DIMENSIONS));
    }
}

void deallocate_bats(struct bat *bats, struct bat *best, struct bat *candidate)
{

    /* for (int i = 0; i < BATS_COUNT; i++) { */
    /*     free(bats[i].position); */
    /*     free(bats[i].velocity); */
    /* } */
    free(bats);

    /* free(best->position); */
    /* free(best->velocity); */
    free(best);

    /* free(candidate->position); */
    /* free(candidate->velocity); */
    free(candidate);
}


void log_bat(struct bat *bat)
{
    logger(LOG_SCALAR_ATRIBUTES, "F,PR,L: %E %E %E\n", bat->frequency, bat->pulse_rate, bat->loudness);

    for (int i = 0; i < DIMENSIONS; i++) {
        logger(LOG_VECTOR_ATRIBUTES, "%E\t%E\t%E\n", bat->velocity[i], bat->position[i], bat->fitness);
    }
}

struct bat get_best(struct bat *bats, struct bat *best)
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
    memcpy(best, &bats[best_indice], sizeof(struct bat));
}

void force_boundry_on_value(double* value)
{
    BOUNDRY_COUNT++;
    if (*value > BOUNDRY_MAX) {
        /* printf("MAX: %E \n", *value); */
        *value = BOUNDRY_MAX;
        BOUNDRY_SCAPE_COUNT++;
        return;
    }
    if (*value < BOUNDRY_MIN) {
        /* printf("MIN: %E\n", *value); */
        *value = BOUNDRY_MIN;
        BOUNDRY_SCAPE_COUNT++;
    }
}

void force_boundry_on_vector(double vector[])
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        force_boundry_on_value(&vector[i]);
    }
}

void update_velocity(struct bat *bat, struct bat *best)
{
    for (int i = 0; i < DIMENSIONS; ++i) {
        bat->velocity[i]+= (bat->position[i] - best->position[i]) * bat->frequency;
        /* printf("Velocity: %E\n", bat->velocity[i]); */
        force_boundry_on_value(&bat->velocity[i]);
    }
}

double generate_frequency()
{
    double beta = my_rand(BETA_MIN, BETA_MAX);
    return FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta;
}

void update_position(struct bat *bat)
{
    for (int i = 0; i < DIMENSIONS; ++i) {
        bat->position[i] += bat->velocity[i];

        force_boundry_on_value(&bat->position[i]);
    }
}

void local_search(struct bat *bat, struct bat *best, double loudness_average)
{
    //double tmp;
    for (int i = 0; i < DIMENSIONS; i++ ) {
        //tmp=best->position[i];
        bat->position[i] = best->position[i] + loudness_average * my_rand(-1.0,1.0);
        /* printf("Position from: %E, Position to %E\n", tmp, bat->position[i]); */
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

void decrease_loudness(struct bat *bat, int iteration)
{
    bat->loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
}

void position_perturbation(struct bat *bat)
{
    int dimension = my_rand(0, DIMENSIONS);
    bat->position[dimension] = bat->position[dimension] * my_rand(0.0,1.0);
    force_boundry_on_vector(bat->position);
}

struct bat* get_worst(struct bat *bats)
{
    double current_worst_val;
    double current_val;

    current_val = current_worst_val = bats[0].fitness;
    int worst_indice = 0;
    for (int i = 0; i < BATS_COUNT; i++) {
        current_val = bats[i].fitness;
        if (current_worst_val <  current_val) {
            current_worst_val = current_val;
            worst_indice = i;
        }
    }

    return &bats[worst_indice];
}

void copy_bat(struct bat *from, struct bat *to)
{
    memcpy(to, from, sizeof(struct bat));
    memcpy(to->position, from->position, sizeof(double) * DIMENSIONS);
    memcpy(to->velocity, from->velocity, sizeof(double) * DIMENSIONS);
}


