#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include "bat.h"

#define ITERATIONS 1000
#define BATS_COUNT 768
#define DIMENSIONS 1000

#define INITIAL_LOUDNESS 1.0

//probability of accepting bad results
#define ALFA 0.5
//affects local search
#define LAMBDA 0.1

#define BETA_MAX 1.0
#define BETA_MIN 0.0

const int EVALUTAION_FUNCTION = GRIEWANK;

const int LOG_OBJECTIVE_ENABLED=1;
const int LOG_ATRIBUTES_ENABLED=1;
const int LOG_RANDOM_ENABLED=0;

void log_bat_stdout(struct bat *bat, int dimensions);
int run_bats(void);
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
void initialize_function(void);

double sphere (double *solution, int dimensions);
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

    initialize_function();

    FREQUENCY_MIN=BOUNDRY_MIN;
    FREQUENCY_MAX=BOUNDRY_MAX;

    my_seed();

    bats = (struct bat *) malloc(sizeof(struct bat) * BATS_COUNT);
    best = (struct bat *) malloc(sizeof(struct bat));
    candidate = (struct bat *) malloc(sizeof(struct bat));
    initialize_bats(bats, best, candidate);

    get_best(bats, best);

    for (iteration = 0; iteration < ITERATIONS ; ++iteration) {
        for (int j = 0; j < BATS_COUNT; ++j){
            bats[j].frequency = generate_frequency();
            update_velocity(&bats[j], best);
            copy_bat(&bats[j], candidate);

            update_position(candidate);

            if (my_rand(0.0,1.0) < candidate->pulse_rate) {
                local_search(candidate, best, calc_loudness_average(bats));
            }

            position_perturbation(candidate);

            bats[j].fitness = fabs(objective_function(bats[j].position, DIMENSIONS));
            candidate->fitness = fabs(objective_function(candidate->position, DIMENSIONS));
            if (my_rand(0.0,1.0) < bats[j].loudness || candidate->fitness < bats[j].fitness) {
                copy_bat(candidate, &bats[j]);
                bats[j].fitness = candidate->fitness;
                bats[j].pulse_rate = 1 - exp(-LAMBDA*iteration);

                decrease_loudness(&bats[j], iteration);
            }
            get_best(bats, best);
            if (LOG_ATRIBUTES_ENABLED) {
                log_bat(&bats[j]);
            }
        }

        if (LOG_OBJECTIVE_ENABLED) {
            average_result = fitness_average(bats);
            worst_result = get_worst(bats)->fitness;
            logger(
                    LOG_OBJECTIVE,
                    "%E\t%E\t%E\n",
                    best->fitness,
                    average_result
                    /* worst_result */
                  );

        }
    }

    log_bat_stdout(best, DIMENSIONS);

    deallocate_bats(bats, best, candidate);
    deallocate_resources();
    return 0;
}


void log_bat_stdout(struct bat *bat, int dimensions) 
{
    logger(LOG_STDOUT, "Best BAT \n");
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        /* logger(LOG_STDOUT, "[%i] = %f\n", i, bat->position[i]); */
        position_average+=bat->position[i];
    }
    position_average/=dimensions;
    /* logger(LOG_STDOUT, "Frequency: %E\n", bat->frequency); */
    /* logger(LOG_STDOUT, "Pulse-rate: %E\n", bat->pulse_rate); */
    /* logger(LOG_STDOUT, "Loudness: %E\n", bat->loudness); */
    /* logger(LOG_STDOUT, "Position Average: %f\n", position_average); */
    printf("ITERATIONS: %d\n", ITERATIONS);
    printf("BATS_COUNT: %d\n", BATS_COUNT);
    printf("DIMENSIONS: %d\n", DIMENSIONS);
    logger(LOG_STDOUT, "Fitness E: %E\n", bat->fitness);
    /* logger(LOG_STDOUT, "Fitness F: %f\n", bat->fitness); */
}

void initialize_bats(struct bat *bats, struct bat *best, struct bat *candidate)
{

    for (int i = 0; i < BATS_COUNT; i++) {
        bats[i].pulse_rate = 0.0;
        bats[i].frequency = 0.0;
        bats[i].fitness;
        bats[i].loudness = INITIAL_LOUDNESS;

        bats[i].velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
        bats[i].position = (double *) malloc(sizeof(double) * DIMENSIONS);
        for (int j = 0; j < DIMENSIONS; j++) {
            bats[i].velocity[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
            bats[i].position[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
        }

        bats[i].fitness = fabs(objective_function(bats[i].position, DIMENSIONS));
    }


    best->velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
    best->position = (double *) malloc(sizeof(double) * DIMENSIONS);

    candidate->velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
    candidate->position = (double *) malloc(sizeof(double) * DIMENSIONS);


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
        logger(LOG_VECTOR_ATRIBUTES, "%E\t%f\t%f\n", bat->velocity[i], bat->position[i], bat->fitness);
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
    copy_bat(&bats[best_indice], best);
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
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = best->position[i] + loudness_average * my_rand(-1.0, 1.0);
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
    // Se a diminuição for lenta o suficiente, o valor encontrado no
    // final da execução tende a ser o ótimo global
    bat->loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);

    /* bat->loudness = INITIAL_LOUDNESS - (INITIAL_LOUDNESS / ITERATIONS) * iteration; */
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
}


void initialize_function(void)
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

