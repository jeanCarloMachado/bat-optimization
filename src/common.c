#include "bat.h"
#include <math.h>
#include <stdarg.h>
#include "mersenne.h"

extern const int LOG_OBJECTIVE_ENABLED;
extern const int LOG_ATRIBUTES_ENABLED;
extern const int LOG_RANDOM_ENABLED;

FILE *LOG_OBJECTIVE_FILE;
FILE *LOG_SCALAR_ATRIBUTES_FILE;
FILE *LOG_VECTOR_ATRIBUTES_FILE;
FILE *LOG_RANDOM_FILE;
int RUN_TIME;

#define DUMP_DIR "./dump"

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

void my_seed(void)
{
    MT_seed();
}

double my_rand(int min, int max)
{

    double result = (double)min + ((max - min)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));

    return result;
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


void log_bat_stdout(struct bat *bat, int dimensions) 
{
    logger(LOG_STDOUT, "Best BAT");
    double position_average =  0;
    for (int i = 0; i < dimensions; i++) {
        logger(LOG_STDOUT, "[%i] = %f\n", i, bat->position[i]);
        position_average+=bat->position[i];
    }
    position_average/=dimensions;
    logger(LOG_STDOUT, "Frequency: %E\n", bat->frequency);
    logger(LOG_STDOUT, "Pulse-rate: %E\n", bat->pulse_rate);
    logger(LOG_STDOUT, "Loudness: %E\n", bat->loudness);
    logger(LOG_STDOUT, "Position Average: %f\n", position_average);
    logger(LOG_STDOUT, "Fitness: %f\n", bat->fitness);
}


