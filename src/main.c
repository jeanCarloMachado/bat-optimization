#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define MAX_ITERATIONS 100
/* #define BATS 40 */
#define BATS_COUNT 40

//pulse rate and frequency are correlated
#define FREQUENCY_MIN 0
#define FREQUENCY_MAX 100

//loudness decrease over time (A)
#define LOUDNESS_MIN 1
#define LOUDNESS_MAX 100


//BAT PARAMETERS
//simillar to the cooling schedule in simulated annealing
#define ALFA 0.9
#define LAMBDA 0.9

#define DIMENSIONS 2

struct bat {
	double pulse_rate; //or frequency
	double loudness;
	double frequency;
	double position[DIMENSIONS];
	double velocity[DIMENSIONS];
};

/* double average_loudness = 0; */
/* bat** solutions; */


/* double current_position() */
/* { */
/*     return previous_position + current_velocity(); */
/* } */

void initialize_bats(struct bat bats[]);
double my_random(double, double);
void my_seed(void);
void print_bat(struct bat bat);
struct bat get_best(struct bat bats[]);
double sphere(double x[], double d);
void print_bat_collection(struct bat bats[]);
double objective_function (struct bat bat);
void update_velocity(struct bat *bat, struct bat best);
double generate_frequency();
void update_position(struct bat *bat);
void local_search(struct bat *bat, struct bat best, double loudness_average);
double calc_loudness_average(struct bat bats[]);

int main() {

	struct bat bats[BATS_COUNT];
	struct bat best;
	struct bat candidate;
	struct bat *current;

	my_seed();

	initialize_bats(bats);

	best = get_best(bats);	

	for (int iteration = 0; iteration < MAX_ITERATIONS ; iteration ++) {
		for (int j = 0; j < BATS_COUNT; j++) {
			current = &bats[j];
			current->frequency = generate_frequency(current);
			update_velocity(current, best);
			candidate = *current;

			update_position(&candidate);

			if (my_random(0,1) < candidate.pulse_rate) {
				local_search(&candidate, best, calc_loudness_average(bats));
			}

			//generate a new solution by flying randomly
			if (my_random(0,1) < bats[j].loudness || objective_function(candidate) < objective_function(*current)) {
				
				memcpy(current->position, candidate.position, sizeof candidate.position);
				current->pulse_rate = 1 - exp(-LAMBDA*iteration);
				current->loudness =  ALFA*current->loudness;
			}

			best = get_best(bats);	

		}
	}



	best = get_best(bats);
	print_bat(best);
	/* //results and vizualization */
	return 0;
}

void local_search(struct bat *bat, struct bat best, double loudness_average)
{
	for (int i = 0; i < DIMENSIONS; i++ ) {
		bat->position[i] = best.position[i] + loudness_average * my_random(0,1);
	}
}

double calc_loudness_average(struct bat bats[])
{
	double total = 0;


	for(int i=0;i<BATS_COUNT;i++) {
		total+= bats[i].loudness;
	}

	return total / BATS_COUNT;
}


void update_velocity(struct bat *bat, struct bat best)
{
	for (int i = 0; i < DIMENSIONS; i++ ) {
		bat->velocity[i] = bat->velocity[i] + (bat->position[i] - best.position[i]) * bat->frequency;
	}
}

void update_position(struct bat *bat)
{
	for (int i = 0; i < DIMENSIONS; i++ ) {
		bat->position[i] = bat->position[i] + bat->velocity[i];
	}
}

double generate_frequency()
{
	double beta = my_random(0,1);
	return (FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN)) * beta;
}


void print_bat_collection(struct bat bats[])
{
	for(int i=0;i<BATS_COUNT;i++) {

		print_bat(bats[i]);
	}

}

void print_bat(struct bat bat)
{
	printf("=== BAT === \n");
	printf("Frequency: %f\n", bat.frequency);
	printf("Loudness: %f\n", bat.loudness);
	printf("Pulse-rate: %f\n", bat.pulse_rate);

	printf("Velocity:\n");
	printf("[0] %f \n", bat.velocity[0]);
	printf("[1] %f \n", bat.velocity[1]);

	printf("Position:\n");
	printf("[0] %f \n", bat.position[0]);
	printf("[1] %f \n", bat.position[1]);
}

struct bat get_best(struct bat bats[])
{
	double current_best_val; 
	double current_val;

 	current_val = current_best_val = objective_function(bats[0]);
	struct bat current_best_bat = bats[0];
	for (int i = 0; i < BATS_COUNT; i++) {
		current_val = objective_function(bats[i]);
		if (current_val < current_best_val) {
			current_best_bat = bats[i];
			current_best_val = current_val;
		}
	}

	return current_best_bat;
}

void my_seed(void)
{
	srand(time(NULL));
}

double my_random(double inferior, double superior)
{
	if (inferior == 0 && superior == 1) {
		return (double)rand() / (double)RAND_MAX ;
	}

	return rand () % (int) superior;
}

void initialize_bats(struct bat bats[])
{
	for (int i = 0; i < BATS_COUNT; i ++ ) {
		bats[i].pulse_rate = 0;
		bats[i].frequency = 0;
		bats[i].loudness = 1;
		bats[i].velocity[0] = 0;
		bats[i].velocity[1] = 0;
		bats[i].position[0] = my_random(0, RAND_MAX);
		bats[i].position[1] = my_random(0, RAND_MAX);
	}
}


double objective_function (struct bat bat)
{
	double result = sphere(bat.position, DIMENSIONS);
	return result;
}

double sphere(double x[], double d)
{
	double total = 0;

	for (int i = 0; i < DIMENSIONS; i++) {
		total+= x[i] * x[i];
	}

	return total;
}

