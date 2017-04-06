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
    my_seed();
    initialize_function(SPHERE);
}

void tearDown()
{
}

void test_initialize_bat(void)
{
    Bat *bat;

    bat = bat_factory();
    TEST_ASSERT(bat_loudness_get(bat) == 1);
    TEST_ASSERT(bat_pulse_rate_get(bat) == 0);
}

void test_copy_bat(void)
{
    Bat *bat_a;
    Bat *bat_b;

    bat_a = bat_factory();
    bat_b = bat_factory();

    TEST_ASSERT(bat_fitness_get(bat_a) != bat_fitness_get(bat_b));
    bat_copy(bat_a, bat_b);

    TEST_ASSERT(bat_fitness_get(bat_a) == bat_fitness_get(bat_b));
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
    RUN_TEST(test_copy_bat);

    return UNITY_END();
}

