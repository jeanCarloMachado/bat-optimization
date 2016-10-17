#include <stdio.h>

#define MAX_ITERATIONS 100
#define BATS 40

#define FREQUENCY_MIN 0
#define FREQUENCY_MAX 100

//loudness decrease over time
//a min
#define LOUDNESS_MIN = 1
//a max
#define LOUDNESS_MAX = 100


//simillar to the cooling schedule in simulated annealing
#define ALFA = 0.9
#define GAMMA = 0.9

int average_loudness = 0;

struct bat {
   float pulse_rate;
   float loudness;
};


float current_frequency()
{
    float beta = random_vector; //0-1
    return FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * BETA
}

float current_velocity()
{
    return  previous_velocity + (previous_position - best_position) * current_frequency()
}


float current_position()
{
    return previous_position + current_velocity();
}

void initialize_bats()
{
    for (int i = 0; i < BATS; i ++ ) {
    //set random loudness and pulse emission
    //loudness tipically 1-2
    //emission rate around 0
    }
}

int main() {


    initialize_bats();
    //  give a random frequency for each
    for (int i = 0; i < MAX_ITERATIONS ; i ++) {
        //generate new solution by adjusting frequency and updating velocities 

        if (rand > ri) {
            //select a solution among the best solutions
            //generate a local solution
        }

        //generate a new solution by flying randomly
        if (rand < Ai && solution > bestSolution) {
            //accept new solutions 
            //increase ri and reduce Ai
        }

        //rank the bats and find the current best
    }
    //results and vizualization
    return 0;
}



