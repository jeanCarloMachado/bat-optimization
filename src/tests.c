#include <stdio.h>
#include "unity.h"

void setUp()
{

}

void tearDown()
{
}

int sum(int a, int b)
{
    return a + b;
}

void test_should_sum()
{
    TEST_ASSERT(sum(2,2) == 4);
}

int main() {

    UNITY_BEGIN();

    RUN_TEST(test_should_sum);

    return UNITY_END();
}

