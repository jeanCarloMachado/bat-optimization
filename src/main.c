#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "mersenne.h"

#define BOUNDRY_MIN -10
#define BOUNDRY_MAX 10

#define DIMENSIONS 256
#define MAX_ITERATIONS 1000
#define BATS_COUNT 40
#define FREQUENCY_MIN 0
#define FREQUENCY_MAX 1
#define LOUDNESS_MIN 0
#define LOUDNESS_MAX 1
#define INITIAL_LOUDNESS 1

#define DUMP_DIR "/home/jean/projects/bat-optimization/dump"
#define BETA_MIN 0
#define BETA_MAX 1

//probability of accepting bad results
#define ALFA 0.1
//affects local search
#define LAMBDA 0.1

#define LOG_OBJECTIVE_ENABLED 1
#define LOG_GENERAL_ENABLED 0
#define LOG_RANDOM_ENABLED 0

#define LOG_OBJECTIVE_FUNCTION 1
#define LOG_FILE_RANDOM 2
#define LOG_STDOUT 3
#define LOG_FILE_GENERAL 4

int RUN_TIME;
FILE *LOG_OBJECTIVE;
FILE *LOG_GENERAL;
FILE *LOG_RANDOM;

struct bat {
    //tends towards 1
    double pulse_rate;
    //tends towards 0
    double loudness;
    double frequency;
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];
};

struct bat get_worst(struct bat bats[]);
void initialize_bats(struct bat bats[]);
double my_rand(int, int);
void my_seed(void);
void log_bat(int destination, struct bat bat);
struct bat get_best(struct bat bats[]);
void update_velocity(struct bat *bat, struct bat best);
double generate_frequency();
void update_position(struct bat *bat);
void local_search(struct bat *bat, struct bat best, double loudness_average);
double calc_loudness_average(struct bat bats[]);
struct bat get_average(struct bat bats[]);
void logger(int destination, char *fmt, ...);
void allocate_resources(void);
void deallocate_resources();
void decrease_loudness(struct bat*, int);
void position_perturbation(struct bat *bat);

double objective_function (struct bat bat);
double sphere(double x[]);
double rastringin (double solution[]);
double griewank (double solution[]);
double ackley (double solution[]);

int main()
{
    allocate_resources();
    struct bat bats[BATS_COUNT];
    struct bat best;
    struct bat worst;
    struct bat average;
    struct bat candidate;
    struct bat *current;
    int iteration;
    double best_result,average_result,worst_result;

    my_seed();

    initialize_bats(bats);
    best = get_best(bats);

    for (iteration = 0; iteration < MAX_ITERATIONS ; iteration ++) {
        for (int j = 0; j < BATS_COUNT; j++) {
            current = &bats[j];
            current->frequency = generate_frequency(current);
            update_velocity(current, best);
            candidate = *current;

            update_position(&candidate);

            if (my_rand(0,1) < candidate.pulse_rate) {
                local_search(&candidate, best, calc_loudness_average(bats));
            }

            position_perturbation(&candidate);

            if (my_rand(0,1) < bats[j].loudness || objective_function(candidate) < objective_function(*current)) {
                memcpy(current->position, candidate.position, sizeof candidate.position);
                current->pulse_rate = 1 - exp(-LAMBDA*iteration);
                decrease_loudness(current, iteration);
            }
            best = get_best(bats);
            if (LOG_GENERAL_ENABLED) {
                log_bat(LOG_FILE_GENERAL, bats[j]);
            }
        }

        best_result = objective_function(best);

        if (LOG_OBJECTIVE_ENABLED) {
            average_result = objective_function(get_average(bats));
            worst_result = objective_function(get_worst(bats));
            logger(
                LOG_OBJECTIVE_FUNCTION,
                "%f\t%f\t%f\n",
                iteration,
                best_result,
                average_result,
                worst_result
            );

        }

        if (fabs(best_result - 0) < 0.0009) {
            break;
        }

    }

    logger(
        LOG_STDOUT,
        "Best of All: %f iterations (%d)",
        objective_function(best),
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
        LOG_OBJECTIVE = fopen(fileName,"w");
        if (LOG_OBJECTIVE == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Objective log: %s\n", fileName);
    }

    if (LOG_GENERAL_ENABLED) {
        sprintf(fileName, "%s/%i-general", DUMP_DIR, RUN_TIME);
        LOG_GENERAL = fopen(fileName,"w");
        if (LOG_GENERAL == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("General log: %s\n", fileName);
    }

    if (LOG_RANDOM_ENABLED) {
        sprintf(fileName, "%s/%i-random", DUMP_DIR, RUN_TIME);
        LOG_RANDOM = fopen(fileName,"w");
        if (LOG_RANDOM == NULL)
        {
            printf("Error opening file %s !\n", fileName);
            exit(1);
        }
        printf ("Random log: %s\n", fileName);
    }
}

void initialize_bats(struct bat bats[])
{
    for (int i = 0; i < BATS_COUNT; i ++ ) {
        bats[i].pulse_rate = 0;
        bats[i].frequency = 0;
        bats[i].loudness = INITIAL_LOUDNESS;

        for (int j = 0; j < DIMENSIONS; j++) {
            bats[i].velocity[j] = 0;
            bats[i].position[j] = my_rand(BOUNDRY_MIN, BOUNDRY_MAX);
        }
    }
}


void deallocate_resources()
{
    if (LOG_OBJECTIVE_ENABLED) {
        fclose(LOG_OBJECTIVE);
    }
    if (LOG_GENERAL_ENABLED) {
        fclose(LOG_GENERAL);
    }
    if (LOG_RANDOM_ENABLED) {
        fclose(LOG_RANDOM);
    }
}

void logger(int destination, char *fmt, ...)
{
    char formatted_string[6666];

    va_list argptr;
    va_start(argptr,fmt);
    vsprintf(formatted_string, fmt, argptr);
    va_end(argptr);

    if (destination == LOG_OBJECTIVE_FUNCTION) 
        fprintf(LOG_OBJECTIVE,"%s",formatted_string);
    else if (destination == LOG_FILE_GENERAL)
        fprintf(LOG_GENERAL,"%s",formatted_string);
    else if (destination == LOG_FILE_RANDOM)
        fprintf(LOG_RANDOM,"%s",formatted_string);
    else if (destination == LOG_STDOUT) 
        printf("%s",formatted_string);
}

void position_perturbation(struct bat *bat)
{
    int dimension = my_rand(0, DIMENSIONS);
    bat->position[dimension] = bat->position[dimension] * my_rand(0,1);
}


void local_search(struct bat *bat, struct bat best, double loudness_average)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = best.position[i] + loudness_average * my_rand(0,1);
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


void update_velocity(struct bat *bat, struct bat best)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->velocity[i] = bat->velocity[i] + (bat->position[i] - best.position[i]) * bat->frequency;
    }
}

void decrease_loudness(struct bat *bat, int iteration)
{
    bat->loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
}

void update_position(struct bat *bat)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = bat->position[i] + bat->velocity[i];
    }
}

double generate_frequency()
{
    double beta = my_rand(BETA_MIN, BETA_MAX);
    return FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta;
}

void log_bat(int destination, struct bat bat)
{
    logger(destination, "F,PR,L: %f %f %f\n", bat.frequency, bat.pulse_rate, bat.loudness);
    /* logger(destination, "\tVelocity:\n"); */
    /* for (int i = 0; i < DIMENSIONS; i++) { */
    /*     logger(destination, "\t[%d] %f \n", i, bat.velocity[i]); */
    /* } */

    /* logger(destination, "\tPosition:\n"); */
    /* for (int i = 0; i < DIMENSIONS; i++) { */
    /*     logger(destination, "\t[%d] %f \n", i, bat.position[i]); */
    /* } */
}


struct bat get_average(struct bat bats[])
{
    struct bat average;

    for (int i = 0; i < BATS_COUNT; i++) {
        average.frequency+=  bats[i].frequency;
        average.loudness+=  bats[i].loudness;
        average.pulse_rate+=  bats[i].pulse_rate;
        for (int j = 0; j < DIMENSIONS; j++) {
            average.position[j]+= bats[i].position[j];
            average.velocity[j]+= bats[i].velocity[j];
        }
    }


    average.frequency =  average.frequency / BATS_COUNT;
    average.loudness =  average.loudness / BATS_COUNT;
    average.pulse_rate =  average.pulse_rate / BATS_COUNT;
    for (int j = 0; j < DIMENSIONS; j++) {
        average.position[j] = average.position[j] / BATS_COUNT;
        average.velocity[j] = average.velocity[j] / BATS_COUNT;
    }

    return average;
}

struct bat get_worst(struct bat bats[])
{
    double current_worst_val;
    double current_val;

    current_val = current_worst_val = objective_function(bats[0]);
    struct bat current_worst_bat = bats[0];
    for (int i = 0; i < BATS_COUNT; i++) {
        current_val = objective_function(bats[i]);
        if (current_val > current_worst_val) {
            current_worst_bat = bats[i];
            current_worst_val = current_val;
        }
    }

    return current_worst_bat;
}



struct bat get_best(struct bat bats[])
{
    double current_best_val; 
    double current_val;

    current_val = current_best_val = objective_function(bats[0]);
    struct bat current_best_bat = bats[0];
    for (int i = 0; i < BATS_COUNT; i++) {
        current_val = objective_function(bats[i]);
        if (current_val < current_best_val) {
            current_best_bat = bats[i];
            current_best_val = current_val;
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
      logger(LOG_FILE_RANDOM, "%i-%i: %f\n", min, max, result); 
  }

   return result;
}

double objective_function (struct bat bat)
{
    double result = sphere(bat.position);
    return result;
}

//best: 0.632
double sphere (double solution[])
{
    double total = 0;

    for (int i = 0; i < DIMENSIONS; i++) {
        total+= solution[i] * solution[i];
    }

    return total;
}

double rastringin (double solution[])
{
    double total = 0;

    for(int i=0;i<DIMENSIONS;i++)
    {
        total=total+(pow(solution[i],(double)2)-10*cos(2*M_PI*solution[i])+10);
    }

   return total;
}

double griewank (double solution[])
{
    double total = 0;

    double top1=0;
    double top2=1;

    for(int i=0;i<DIMENSIONS;i++)
    {
        top1=top1+pow((solution[i]),(double)2);
        top2=top2*cos((((solution[i])/sqrt((double)(i+1)))*M_PI)/180);
    }
    total=(1/(double)4000)*top1-top2+1;

   return total;
}

double ackley (double solution[])
{
    int i;
    double aux1, aux;

    for (i = 0; i < DIMENSIONS; i++)
    {
        aux += solution[i]*solution[i];
    }
    for (i = 0; i < DIMENSIONS; i++)
    {
        aux1 += cos(2.0*M_PI*solution[i]);
    }

    return (-20.0*(exp(-0.2*sqrt(1.0/(double)DIMENSIONS*aux)))-exp(1.0/(double)DIMENSIONS*aux1)+20.0+exp(1));
}
