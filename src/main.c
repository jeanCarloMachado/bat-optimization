#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#define DIMENSIONS 2
#define MAX_ITERATIONS 100000
#define BATS_COUNT 40
#define FREQUENCY_MIN 0
#define FREQUENCY_MAX 100
#define LOUDNESS_MIN 1
#define LOUDNESS_MAX 100
#define ALFA 0.5
#define LAMBDA 0.1
#define DUMP_DIR "/home/jean/projects/bat-optimization/dump"

#define DEBUG_LEVEL 1
#define DEBUG_RANDOM 0

#define LOG_FILE_MAIN 1
#define LOG_FILE_RANDOM 2
#define LOG_STDOUT 3

int RUN_TIME;
FILE *LOG;
FILE *LOG_RANDOM;

struct bat {
    double pulse_rate;
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
double sphere(double x[], double d);
double objective_function (struct bat bat);
void update_velocity(struct bat *bat, struct bat best);
double generate_frequency();
void update_position(struct bat *bat);
void local_search(struct bat *bat, struct bat best, double loudness_average);
double calc_loudness_average(struct bat bats[]);
struct bat get_average(struct bat bats[]);
void logger(int destination, char *fmt, ...);
void allocate_resources(void);
void deallocate_resources();

int main()
{
    allocate_resources();
    struct bat bats[BATS_COUNT];
    struct bat best;
    struct bat worst;
    struct bat average;
    struct bat candidate;
    struct bat *current;

    my_seed();

    initialize_bats(bats);
    best = get_best(bats);

    for (int iteration = 0; iteration < MAX_ITERATIONS ; iteration ++) {
        logger(LOG_FILE_MAIN, "Iteration %i\n", iteration);
        for (int j = 0; j < BATS_COUNT; j++) {
            current = &bats[j];
            current->frequency = generate_frequency(current);
            update_velocity(current, best);
            candidate = *current;

            update_position(&candidate);

            if (my_rand(0,1) < candidate.pulse_rate) {
                local_search(&candidate, best, calc_loudness_average(bats));
                if (DEBUG_LEVEL >= 2) {
                    logger(LOG_FILE_MAIN, "Doing Local Search\n");
                }
            }

            if (my_rand(0,1) < bats[j].loudness || objective_function(candidate) < objective_function(*current)) {
                memcpy(current->position, candidate.position, sizeof candidate.position);
                current->pulse_rate = 1 - exp(-LAMBDA*iteration);
                current->loudness =  ALFA*current->loudness;
                if (DEBUG_LEVEL >= 2) {
                    logger(LOG_FILE_MAIN, "Updating with local search\n");
                }
            }
            best = get_best(bats);
            worst = get_worst(bats);
        }

        logger(LOG_FILE_MAIN, "Iteration Best\n");
        log_bat(LOG_FILE_MAIN, best);
        logger(LOG_FILE_MAIN, "Iteration Average\n");
        average = get_average(bats);
        log_bat(LOG_FILE_MAIN, average);
    }

    logger(LOG_STDOUT, "BEST");
    log_bat(LOG_STDOUT, best);
    logger(LOG_STDOUT, "AVERAGE");
    average = get_average(bats);
    log_bat(LOG_STDOUT, average);
    logger(LOG_STDOUT, "WORST");
    log_bat(LOG_STDOUT, worst);

    deallocate_resources();
    return 0;
}

void allocate_resources()
{
    RUN_TIME = time(NULL);

    char fileName[100];
    sprintf(fileName, "%s/%i-main", DUMP_DIR, RUN_TIME);
    LOG = fopen(fileName,"w");
    if (LOG == NULL)
    {
        printf("Error opening file %s !\n", fileName);
        exit(1);
    }
    printf ("Main log: %s\n", fileName);

    if (DEBUG_RANDOM) {
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

void deallocate_resources()
{
    fclose(LOG);
    if (DEBUG_RANDOM) {
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

    /* printf("%s",formatted_string); */

    if (destination == LOG_FILE_MAIN) 
        fprintf(LOG,"%s",formatted_string);
    else if (destination == LOG_FILE_RANDOM)
        fprintf(LOG_RANDOM,"%s",formatted_string);
    else if (destination == LOG_STDOUT) 
        printf("%s",formatted_string);
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

void update_position(struct bat *bat)
{
    for (int i = 0; i < DIMENSIONS; i++ ) {
        bat->position[i] = bat->position[i] + bat->velocity[i];
    }
}

double generate_frequency()
{
    double beta = my_rand(0,1);
    return (FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN)) * beta;
}

void log_bat(int destination, struct bat bat)
{
    logger(destination, " = BAT = \n");
    logger(destination, "\tFrequency: %f\n", bat.frequency);
    logger(destination, "\tLoudness: %f\n", bat.loudness);
    logger(destination, "\tPulse-rate: %f\n", bat.pulse_rate);

    logger(destination, "\tVelocity:\n");
    for (int i = 0; i < DIMENSIONS; i++) {
        logger(destination, "\t[%d] %f \n", i, bat.velocity[i]);
    }

    logger(destination, "\tPosition:\n");
    for (int i = 0; i < DIMENSIONS; i++) {
        logger(destination, "\t[%d] %f \n", i, bat.position[i]);
    }
}


struct bat get_average(struct bat bats[])
{
    struct bat average;

    for (int i = 0; i < BATS_COUNT; i++) {
        average.frequency =  bats[i].frequency;
        average.loudness =  bats[i].loudness;
        average.pulse_rate =  bats[i].pulse_rate;
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
    srand(time(NULL));
}

double my_rand(int min, int max)
{
    double result;

    if (min == 0 && max == 1) {
        result = (double)rand() / (double)RAND_MAX ;

        if (DEBUG_RANDOM) {
            logger(LOG_FILE_RANDOM, "0-1: %f\n", result);
        }
        return result;
    }

    result =  (rand() % (max + 1 - min)) + min;
    if (DEBUG_RANDOM) {
        logger(LOG_FILE_RANDOM, "0-100: %i\n", result);
    }
    return result;
}

void initialize_bats(struct bat bats[])
{
    for (int i = 0; i < BATS_COUNT; i ++ ) {
        bats[i].pulse_rate = 0;
        bats[i].frequency = 0;
        bats[i].loudness = 1;

        for (int j = 0; j < DIMENSIONS; j++) {
            bats[i].velocity[j] = 0;
            bats[i].position[j] = my_rand(0, 100);
        }
    }
}

double objective_function (struct bat bat)
{
    double result = sphere(bat.position, DIMENSIONS);
    return result;
}

double sphere(double x[], double d)
{
    double total = 0;

    for (int i = 0; i < DIMENSIONS; i++) {
        total+= x[i] * x[i];
    }

    return total;
}

