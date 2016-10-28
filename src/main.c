#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#define DIMENSIONS 200
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
#define ALFA 0.9
//affects local search
#define LAMBDA 0.9


#define LOG_LEVEL 1

#define LOG_FILE_MAIN 1
#define LOG_FILE_RANDOM 2
#define LOG_STDOUT 3

int RUN_TIME;
FILE *LOG;
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
        }

        if (LOG_LEVEL >= 1) {
            int best_result = abs(objective_function(best));
            int average_result = abs(objective_function(get_average(bats)));
            logger(
                LOG_FILE_MAIN,
                "%u\t%d\t%u\n",
                iteration,
                best_result,
                average_result
            );

        }
    }

    printf("Best of All: %f", objective_function(best));

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

    if (LOG_LEVEL > 9) {
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
            bats[i].position[j] = my_rand(0, 100);
        }
    }
}


void deallocate_resources()
{
    fclose(LOG);
    if (LOG_LEVEL > 9) {
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

    if (destination == LOG_FILE_MAIN) 
        fprintf(LOG,"%s",formatted_string);
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

        if (LOG_LEVEL > 9) {
            logger(LOG_FILE_RANDOM, "0-1: %f\n", result);
        }
        return result;
    }


    if (min == -1 && max == 1) {
        double signal;
        result = (double)rand() / (double)RAND_MAX ;

        if (LOG_LEVEL > 9) {
            logger(LOG_FILE_RANDOM, "0-1: %f\n", result);
        }

       signal = (double)rand() / (double)RAND_MAX ;
       if (signal < 0.5)  {
            result = -result;
       }

        return result;
    }


    result =  (rand() % (max + 1 - min)) + min;
    if (LOG_LEVEL > 9) {
        logger(LOG_FILE_RANDOM, "0-100: %i\n", result);
    }
    return result;
}

double objective_function (struct bat bat)
{
    double result = griewank(bat.position);
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
