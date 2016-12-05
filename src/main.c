#include <stdio.h>
#include "bat.h"
#include <time.h>

extern int run_bats(void);
int main(void) {
    clock_t begin = clock();
    run_bats();
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time took: %f\n", time_spent);
    return 0;
}

