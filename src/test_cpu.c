#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include "unity.h"
#include "bat/internal.h"


int bats_count;
int evaluation_function;
int iterations;

void setUp()
{
}

void tearDown()
{
}

void test_initialize_bat(void)
{
    Bat *bat;

    bat = bat_factory();
    TEST_ASSERT(bat->loudness == 1);
    TEST_ASSERT(bat->pulse_rate == 0);
}

void test_schewefel(void)
{
    double solution[3];
    solution[0]=499.453719;
    solution[1]=-499.453715;
    solution[2]=-0.000249;

    TEST_ASSERT(fabs(0.000016 - schwefel(solution, 3)) < 000.1);
}

int main(void) {

    UNITY_BEGIN();

    RUN_TEST(test_schewefel);
    RUN_TEST(test_initialize_bat);

    return UNITY_END();
}

