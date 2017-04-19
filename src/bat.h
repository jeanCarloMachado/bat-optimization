#ifndef BAT_H_
#define BAT_H_
#include "bat/internal.h"

typedef void Bat;
extern int bats_count;
extern int evaluation_function;
extern int iterations;

void bat_run(void);

/* struct Bat* bat_factory(); */
/* double bat_loudness_get(struct Bat *bat); */
/* double bat_pulse_rate_get(struct Bat *bat); */
/* double bat_fitness_get(struct Bat *bat); */
/* void bat_copy(struct Bat *from, struct Bat *to); */

#endif
