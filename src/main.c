#include <stdio.h>
#include "bat.h"
#include <time.h>
#include <string.h>

extern int run_bats(void);

int main(int argc, char **argv)
{
    char *HELP = "--help";

    if (argc > 1 && strcmp(argv[1], HELP) == 0) {
        printf("The CPU version of the BAT algorithm");
        return 0;
    }

    clock_t begin = clock();
    run_bats();
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time took: %f\n", time_spent);
    return 0;
}

