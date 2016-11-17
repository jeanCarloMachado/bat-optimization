#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include "unity.h"
#include "bat.h"

int BATS_COUNT = 2;
struct bat *bats;

extern struct bat get_best(struct bat *bats, struct bat *best);

void setUp()
{
    bats = (struct bat *) malloc(sizeof(struct bat) * BATS_COUNT);
}

void tearDown()
{
    free(bats);
}

int sum(int a, int b)
{
    return a + b;
}

void test_get_best()
{
    extern struct bat *bats;
    struct bat *best;

    best = (struct bat *) malloc(sizeof(struct bat));

    for (int i = 0; i < DIMENSIONS; i++ ) {
        bats[i].fitness = 1/(i+1);
    }

    get_best(bats, best);

    TEST_ASSERT(best->fitness == bats[DIMENSIONS-1].fitness);
    free(best);
}

int main() {

    UNITY_BEGIN();

    RUN_TEST(test_get_best);

    return UNITY_END();
}

