#ifndef BAT_H_
#define BAT_H_

#define DIMENSIONS 10

struct bat {
    //tends towards 1
    double pulse_rate;
    //tends towards 0
    double loudness;
    double fitness;
    double frequency;
    double position[DIMENSIONS];
    double velocity[DIMENSIONS];
};

int run_bats();

struct bat get_best(struct bat *bats, struct bat *best);

#endif
