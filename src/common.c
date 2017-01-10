#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include "bat.h"

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

