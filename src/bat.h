#ifndef BAT_H_
#define BAT_H_
#include "bat/internal.h"

struct Bat* bat_factory();
int bats_count;
int evaluation_function;
int iterations;

typedef void Bat;
int bat_run(void);

double bat_loudness_get(struct Bat *bat);
double bat_pulse_rate_get(struct Bat *bat);
double bat_fitness_get(struct Bat *bat);

void bat_copy(struct Bat *from, struct Bat *to);

#endif
