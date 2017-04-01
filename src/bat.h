#ifndef BAT_H_
#define BAT_H_
#include "bat/internal.h"

int bats_count;
int evaluation_function;
int iterations;

typedef void Bat;
int bats_run(void);

#endif
