#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include "unity.h"
#include "bat.h"

void setUp()
{
}

void tearDown()
{
}


void test_copy_bat()
{
    struct bat from;
    struct bat to;
    double position_from[3] = {1.0,2.0,3.0};
    double position_to[3] = {4.0,5.0,6.0};

    from.position = position_from;
    TEST_ASSERT(2.0 == from.position[1]);

    to.position = position_to;
    TEST_ASSERT(5.0 == to.position[1]);

    copy_bat(&from, &to);

    TEST_ASSERT(2.0 == to.position[1]);
}



void test_get_worst()
{
    struct bat *bats;
    struct bat *worst;
    int BATS_COUNT=10;

    bats = (struct bat *) malloc(sizeof(struct bat) * BATS_COUNT);

    for (int i = 0; i < BATS_COUNT; i++ ) {
        bats[i].fitness = 1.0/(i+1);
    }

    worst = get_worst(bats);

    TEST_ASSERT(worst->fitness == bats[0].fitness);
    free(bats);
}

void test_get_best()
{
    struct bat *bats;
    struct bat *best;
    int BATS_COUNT=5;

    bats = (struct bat *) malloc(sizeof(struct bat) * BATS_COUNT);
    best = (struct bat *) malloc(sizeof(struct bat));

    for (int i = 0; i < BATS_COUNT; i++ ) {
        bats[i].fitness = 1/(i+1);
    }

    get_best(bats, best);

    TEST_ASSERT(best->fitness == bats[BATS_COUNT-1].fitness);

    free(best);
    free(bats);
}

void test_rosenbrook()
{
    double solution[3];
    solution[0]=1;
    solution[1]=1;
    solution[2]=1;

    TEST_ASSERT(0.0 == rosenbrock(solution, 3));


    solution[0]=1;
    solution[1]=0;
    solution[2]=0;

    TEST_ASSERT(101 == (int) rosenbrock(solution, 3));
}

int main() {

    UNITY_BEGIN();

    RUN_TEST(test_get_best);
    RUN_TEST(test_get_worst);
    RUN_TEST(test_rosenbrook);
    RUN_TEST(test_copy_bat);

    return UNITY_END();
}

