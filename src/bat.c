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
#define ITERATIONS 10000
#define BATS_COUNT 256
#define DIMENSIONS 100

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

const int EVALUTAION_FUNCTION = ACKLEY;

const int LOG_OBJECTIVE_ENABLED=0;
const int LOG_ATRIBUTES_ENABLED=0;

int BOUNDRY_MAX;
int BOUNDRY_MIN;
int FREQUENCY_MIN;
int FREQUENCY_MAX;
int FIRST_SEED;
int SECOND_SEED;

int BOUNDRY_SCAPE_COUNT = 0;
int BOUNDRY_COUNT = 0;


double (*objective_function)(double[], int);
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

FILE *LOG_OBJECTIVE_FILE;
FILE *LOG_SCALAR_ATRIBUTES_FILE;
FILE *LOG_VECTOR_ATRIBUTES_FILE;
int RUN_TIME;

double shuber (double solution[], int dimensions)
{
    /*
       -   Domain  |x| <= 10.0
       -   Number of local minimum = 400
       -   Global minimum fmin = -24.062499 at the ff. points
       -    (-6.774576, -6.774576), ..., (5.791794, 5.791794)
       */
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
    for (int i = 0; i < BATS_COUNT; i++) {
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
    else if (destination == LOG_SCALAR_ATRIBUTES)
        fprintf(LOG_SCALAR_ATRIBUTES_FILE,"%s",formatted_string);
    else if (destination == LOG_VECTOR_ATRIBUTES)
        fprintf(LOG_VECTOR_ATRIBUTES_FILE,"%s",formatted_string);
    else if (destination == LOG_STDOUT) 
        printf("%s",formatted_string);
}

void my_seed(void)
{
    FIRST_SEED = time(NULL);
    SECOND_SEED = 19;
}

uint64_t xorshift128plus(void) {

    uint64_t x = FIRST_SEED;
    uint64_t const y = SECOND_SEED;
    FIRST_SEED = y;
    x ^= x << 23; // a
    SECOND_SEED = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
    return SECOND_SEED + y;
}

double my_rand(int min, int max)
{
    double result = (double)min + ((max - min)*(double)xorshift128plus()/(RAND_MAX+1.0));
    return result;
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
}

void log_bat_stdout(struct bat *bat, int dimensions) 
{
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        position_average+=bat->position[i];
    }
    position_average/=dimensions;
    printf("ITERATIONS: %d\n", ITERATIONS);
    printf("BATS_COUNT: %d\n", BATS_COUNT);
    printf("DIMENSIONS: %d\n", DIMENSIONS);
    logger(LOG_STDOUT, "Fitness E: %E\n", bat->fitness);
}

void initialize_bats(struct bat *bats, struct bat *best, struct bat *candidate)
{

    for (int i = 0; i < BATS_COUNT; i++) {
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
    free(bats);
    free(best);
    free(candidate);
}


void log_bat(struct bat *bat)
{
    logger(LOG_SCALAR_ATRIBUTES, "F,PR,L: %E %E %E\n", bat->frequency, bat->pulse_rate, bat->loudness);

    for (int i = 0; i < DIMENSIONS; i++) {
        logger(LOG_VECTOR_ATRIBUTES, "%E\t%f\t%f\n", bat->velocity[i], bat->position[i], bat->fitness);
    }
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
    for (int i = 0; i < BATS_COUNT; i++) {
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
        }

    }

    printf("Function %s\n", get_function_name(EVALUTAION_FUNCTION));
    log_bat_stdout(best, DIMENSIONS);

    deallocate_bats(bats, best, candidate);
    deallocate_resources();
    return 0;
}
