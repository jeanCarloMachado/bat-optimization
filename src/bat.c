#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include "bat.h"

#define ALFA 0.5
#define LAMBDA 0.1
#define BETA_MAX 1.0
#define BETA_MIN 0.0
#define INITIAL_LOUDNESS 1.0
#define DIMENSIONS 1000

const int LOG_OBJECTIVE_ENABLED=1;

extern int bats_count;
extern int iterations;
extern int evaluation_function;

struct bat {
    //tends towards 1
    double pulse_rate;
    //tends towards 0
    double loudness;
    double fitness;
    double frequency;
    double *position;
    double *velocity;
};


int boundry_max;
int BOUNDRY_MIN;
int FREQUENCY_MIN;
int FREQUENCY_MAX;
int FIRST_SEED;
int SECOND_SEED;

int BOUNDRY_SCAPE_COUNT = 0;
int BOUNDRY_COUNT = 0;

double (*objective_function)(double[], int);
void bat_stdout(struct bat *bat, int dimensions);
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
void log_objective(struct bat *best, struct bat *all_bats);
void initialize_bats(struct bat *bats, struct bat *best, struct bat *candidate);
void deallocate_bats(struct bat *bats, struct bat *best, struct bat *candidate);
void initialize_function(void);

FILE *LOG_OBJECTIVE_FILE;
FILE *LOG_SCALAR_ATRIBUTES_FILE;
FILE *LOG_VECTOR_ATRIBUTES_FILE;
int RUN_TIME;

double shuber (double solution[], int dimensions)
{
    double sum = 0.0;
    for (int i = 0; i < dimensions; i++) {
        sum += -sin(2.0*solution[i]+1.0)
            -2.0*sin(3.0*solution[i]+2.0)
            -3.0*sin(4.0*solution[i]+3.0)
            -4.0*sin(5.0*solution[i]+4.0)
            -5.0*sin(6.0*solution[i]+5.0);
    }
}


double sphere (double *solution, int dimensions)
{
    double total = 0;

    for (int i = 0; i < dimensions; i++) {
        total+= solution[i] * solution[i];
    }

    return total;
}


double rosenbrock (double solution[], int dimensions)
{
    double total = 0;
    for (int i = 0; i < dimensions-1; i++)
    {
        total=total+100.*pow((solution[i+1] - pow(solution[i],2.)),2) + pow((1. - solution[i]),2);
    }

    return total;
}


double schwefel(double solution[], int dimensions)
{
    double aux = 0;
    for (int i=0;i<dimensions;i++)
    {
        aux += solution[i]*sin(sqrt(fabs(solution[i])));
    }
    return(-1*aux/dimensions);
}

double rastringin (double solution[], int dimensions)
{
    double total = 0;

    for(int i=0;i<dimensions;i++)
    {
        total=total+(pow(solution[i],(double)2)-10*cos(2*M_PI*solution[i])+10);
    }

    return total;
}

double griewank(double solution[], int dimensions)
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

double ackley(double solution[], int dimensions)
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

    result = -20.0*(exp(-0.2*sqrt(1.0/(float)dimensions*aux)))-exp(1.0/(float)dimensions*aux1)+20.0+exp(1);

    return result;
}

void copy_bat(struct bat *from, struct bat *to)
{
    memcpy(to, from, sizeof(struct bat));
}

struct bat get_best(struct bat *bats, struct bat *best)
{
    double current_best_val; 
    int best_indice;

    current_best_val = bats[0].fitness;
    best_indice = 0;
    for (int i = 0; i < bats_count; i++) {
        if (bats[i].fitness < current_best_val) {
            current_best_val = bats[i].fitness;
            best_indice = i;
        }
    }
    copy_bat(&bats[best_indice], best);
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
    else if (destination == LOG_STDOUT)
        printf("%s",formatted_string);
}

void my_seed(void)
{
    /* FIRST_SEED = time(NULL); */
    /* SECOND_SEED = 19; */
    MT_seed();
}

double my_rand( double inferior, double superior)
{
  double result = (double)inferior + ((superior - inferior)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));

    return result;
}


void deallocate_resources()
{
    if (LOG_OBJECTIVE_ENABLED) {
        fclose(LOG_OBJECTIVE_FILE);
    }
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
}

void bat_stdout(struct bat *bat, int dimensions) 
{
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        position_average+=bat->position[i];
    }
    position_average/=dimensions;
    printf("ITERATIONS: %d\n", iterations);
    printf("BATS_COUNT: %d\n", bats_count);
    printf("DIMENSIONS: %d\n", DIMENSIONS);
    printf("POPULATION: %d\n", bats_count);
    logger(LOG_STDOUT, "Fitness E: %E\n", bat->fitness);
}

void initialize_bats(struct bat *bats, struct bat *best, struct bat *candidate)
{

    for (int i = 0; i < bats_count; i++) {
        bats[i].fitness;
        bats[i].loudness = INITIAL_LOUDNESS;

        bats[i].velocity = (double *) malloc(sizeof(double) * DIMENSIONS);
        bats[i].position = (double *) malloc(sizeof(double) * DIMENSIONS);
        for (int j = 0; j < DIMENSIONS; j++) {
            bats[i].velocity[j] = my_rand(BOUNDRY_MIN, boundry_max);
            bats[i].position[j] = my_rand(BOUNDRY_MIN, boundry_max);
            //printf("%f %d\n",my_rand(BOUNDRY_MIN, boundry_max), boundry_max);
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
    free(bats);
    free(best);
    free(candidate);
}


void log_objective(struct bat *best, struct bat *all_bats)
{
    double average_fitness = fitness_average(all_bats);
    logger(LOG_OBJECTIVE, "%E %E\n", best->fitness, average_fitness);
}


void force_boundry_on_value(double* value)
{
    BOUNDRY_COUNT++;
    if (*value > boundry_max) {
        *value = boundry_max;
        BOUNDRY_SCAPE_COUNT++;
        return;
    }
    if (*value < BOUNDRY_MIN) {
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


    for(int i=0;i<bats_count;i++) {
        total+= bats[i].loudness;
    }

    return total / bats_count;
}

double fitness_average(struct bat bats[])
{
    double result = 0;

    for (int i = 0; i < bats_count; i++) {
        result+= bats[i].fitness;
    }

    return result / bats_count;
}

void decrease_loudness(struct bat *bat, int iteration)
{
    bat->loudness = INITIAL_LOUDNESS*pow(ALFA, iteration);
}

void position_perturbation(struct bat *bat)
{
    int dimension = my_rand(0, DIMENSIONS-1);
    //bat->position[dimension] = bat->position[dimension] * my_rand(0.0,1.0);
    force_boundry_on_vector(bat->position); 
}

struct bat* get_worst(struct bat *bats)
{
    double current_worst_val;
    double current_val;

    current_val = current_worst_val = bats[0].fitness;
    int worst_indice = 0;
    for (int i = 0; i < bats_count; i++) {
        current_val = bats[i].fitness;
        if (current_worst_val <  current_val) {
            current_worst_val = current_val;
            worst_indice = i;
        }
    }

    return &bats[worst_indice];
}


void initialize_function(void)
{
    switch(evaluation_function) {
        case SPHERE:
            BOUNDRY_MIN = 0.0;
            boundry_max = 100.0;
            objective_function = &sphere;
            break;
        case RASTRINGIN:
            BOUNDRY_MIN = -5.12;
            boundry_max = 5.12;
            objective_function = &rastringin;
            break;
        case GRIEWANK:
            BOUNDRY_MIN = -600.0;
            boundry_max = 600.0;
            objective_function = &griewank;
            break;
        case ACKLEY:
            BOUNDRY_MIN = -32.0;
            boundry_max = 32.0;
            objective_function = &ackley;
            break;
        case SHUBER:
            BOUNDRY_MIN = -100.0;
            boundry_max = 100.0;
            objective_function = &shuber;
            break;
        case SCHWEFEL:
            BOUNDRY_MIN = -500.0;
            boundry_max = 500.0;
            objective_function = &schwefel;
            break;
        case ROSENBROOK:
            BOUNDRY_MIN = -30.0;
            boundry_max = 30.0;
            objective_function = &rosenbrock;
            break;
    }
}

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
    FREQUENCY_MAX=boundry_max;

    my_seed();

    bats = (struct bat *) malloc(sizeof(struct bat) * bats_count);
    best = (struct bat *) malloc(sizeof(struct bat));
    candidate = (struct bat *) malloc(sizeof(struct bat));
    initialize_bats(bats, best, candidate);

    get_best(bats, best);

    for (iteration = 0; iteration < iterations ; ++iteration) {
        for (int j = 0; j < bats_count; ++j) {
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
            if (my_rand(0.0,1.0) < bats[j].loudness && candidate->fitness < bats[j].fitness) {
                copy_bat(candidate, &bats[j]);
                bats[j].fitness = candidate->fitness;
                bats[j].pulse_rate = 1 - exp(-LAMBDA*iteration);

                decrease_loudness(&bats[j], iteration);
            }
        }
        get_best(bats, best);
        if (LOG_OBJECTIVE_ENABLED) {
            log_objective(best, bats);
        }
    }

    bat_stdout(best,DIMENSIONS);
    printf("Function %s\n", get_function_name(evaluation_function));

    deallocate_bats(bats, best, candidate);
    deallocate_resources();
    return 0;
}
