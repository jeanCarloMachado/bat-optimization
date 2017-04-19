#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "bat.h"

extern void bat_run(void);
int iterations = 100;
int bats_count = 768;
int evaluation_function = 0;

char* get_function_name(int index)
{
    char* str;
    switch (index) {
        case ROSENBROOK:
            str = "ROSENBROOK";
            break;
        case SPHERE:
            str = "SPHERE";
            break;
        case SCHWEFEL:
            str = "SCHWEFEL";
            break;
        case ACKLEY:
            str = "ACKLEY";
            break;
        case RASTRINGIN:
            str = "RASTRINGIN";
            break;
        case GRIEWANK:
            str = "GRIEWANK";
            break;
        case SHUBER:
            str = "SHUBER";
            break;
    }

    return str;
}



int main(int argc, char **argv)
{
    char *HELP = "--help";

    if (argc > 1 && strcmp(argv[1], HELP) == 0) {
        printf("The CPU version of the BAT algorithm\
                You may optionally pass the given variables:\
                ITERATIONS=1000\
                BATS_COUNT=1000\
                EVALUATION_FUNCTION=1\
\
                where EVALUATION_FUNCTION can be one of the following:\
\
                0 ROSENBROOK,\
                1 SPHERE,\
                2 SCHWEFEL,\
                3 ACKLEY,\
                4 RASTRINGIN,\
                5 GRIEWANK,\
                6 SHUBER\
\
                ");
        return 0;
    }

    char* sBatsCount;
    sBatsCount = getenv("BATS_COUNT");
    if (sBatsCount != NULL) {
        bats_count = atoi(sBatsCount);
    }

    char* sIterations;
    sIterations = getenv("ITERATIONS");
    if (sIterations != NULL) {
        iterations = atoi(sIterations);
    }

    char* sEvaluationFunction;
    sEvaluationFunction = getenv("EVALUATION_FUNCTION");
    if (sEvaluationFunction != NULL) {
        evaluation_function = atoi(sEvaluationFunction);
    }

    clock_t begin = clock();
    bat_run();
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time took: %f\n", time_spent);
    printf("Function %s\n", get_function_name(evaluation_function));
    return 0;
}

