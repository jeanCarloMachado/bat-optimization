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

int main() {

    UNITY_BEGIN();

    RUN_TEST(test_get_best);
    RUN_TEST(test_get_worst);

    return UNITY_END();
}

