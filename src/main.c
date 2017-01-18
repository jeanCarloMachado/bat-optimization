#include <stdio.h>
#include <stdlib.h>
#include "bat.h"
#include <time.h>
#include <string.h>

extern int run_bats(void);
int iterations = 10000;
int bats_count = 768;
int evaluation_function = ACKLEY;

int main(int argc, char **argv)
{
    char *HELP = "--help";

    if (argc > 1 && strcmp(argv[1], HELP) == 0) {
        printf("The CPU version of the BAT algorithm\
                You may optionally pass the given variables:\
                ITERATIONS=1000\
                BATS_COUNT=1000\
                FUNCTION_NUM=1\
\
                where FUNCTION_NUM can be one of the following:\
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
    run_bats();
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time took: %f\n", time_spent);
    return 0;
}

