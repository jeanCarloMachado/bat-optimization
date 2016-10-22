#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_ITERATIONS 100
/* #define BATS 40 */
#define BATS 3

//pulse rate and frequency are correlated
#define FREQUENCY_MIN 0
#define FREQUENCY_MAX 100

//loudness decrease over time (A)
#define LOUDNESS_MIN 1
#define LOUDNESS_MAX 100


//BAT PARAMETERS
//simillar to the cooling schedule in simulated annealing
#define ALFA 0.9
#define GAMMA 0.9

#define DIMENSIONS 2
#define SEARCH_BOUNDRY 100

struct bat {
	int pulse_rate; //or frequency
	int loudness;
	int frequency;
	int position[DIMENSIONS];
	int velocity[DIMENSIONS];
};

/* int average_loudness = 0; */
/* bat** solutions; */

/* int current_frequency() */
/* { */
/*     int beta = random_vector; //0-1 */
/*     return FREQUENCY_MIN + (FREQUENCY_MAX - FREQUENCY_MIN) * beta */
/* } */

/* int current_velocity() */
/* { */
/*     return  previous_velocity + (previous_position - best_position) * current_frequency() */
/* } */


/* int current_position() */
/* { */
/*     return previous_position + current_velocity(); */
/* } */

void initialize_bats(struct bat bats[]);
int my_random(int, int);
int my_seed(void);
void print_bat(struct bat bat);
struct bat get_best(struct bat bats[]);
int sphere(int x[], int d);
void print_bat_collection(struct bat bats[]);
double objective_function (struct bat bat);

int main() {

	struct bat bats[BATS];

	my_seed();

	initialize_bats(bats);

	struct bat current_best = get_best(bats);	
	print_bat_collection(bats);
	printf("BEST BAT\n");
	print_bat(current_best);

	/* for (int i = 0; i < MAX_ITERATIONS ; i ++) { */
	/*     //generate new solution by adjusting frequency and updating velocities */ 

	/*     if (rand > ri) { */
	/*         //select a solution among the best solutions */
	/*         //generate a local solution */
	/*     } */

	/*     //generate a new solution by flying randomly */
	/*     if (rand < average_loudness && solution > bestSolution) { */
	/*         //accept new solutions */ 
	/*         //increase ri and reduce Ai */
	/*     } */

	/*     //rank the bats and find the current best */
	/* } */
	/* //results and vizualization */
	return 0;
}


void print_bat_collection(struct bat bats[])
{
	for(int i=0;i<BATS;i++) {

		print_bat(bats[i]);
	}

}

void print_bat(struct bat bat)
{
	printf("Bat Position\n");
	printf("[0]: %i \n", bat.position[0]);
	printf("[1]: %i \n", bat.position[1]);
}

struct bat get_best(struct bat bats[])
{
	int current_best_val; 
	int current_val;

 	current_val = current_best_val = objective_function(bats[0]);
	struct bat current_best_bat = bats[0];
	for (int i = 0; i < BATS; i++) {
		current_val = objective_function(bats[i]);
		if (current_val < current_best_val) {
			current_best_bat = bats[i];
			current_best_val = current_val;
		}
	}

	return current_best_bat;
}

int my_seed(void)
{
	srand(time(NULL));
}

int my_random(int inferior, int superior)
{
	return rand () % superior;
}

void initialize_bats(struct bat bats[])
{
	for (int i = 0; i < BATS; i ++ ) {
		bats[i].pulse_rate = 0;
		bats[i].frequency = 0;
		bats[i].loudness = 1;
		bats[i].velocity[0] = 0;
		bats[i].velocity[1] = 0;
		bats[i].position[0] = my_random(0, SEARCH_BOUNDRY);
		bats[i].position[1] = my_random(0, SEARCH_BOUNDRY);
	}
}

double objective_function (struct bat bat)
{
	int result = sphere(bat.position, DIMENSIONS);
	return result;
}

int sphere(int x[], int d)
{
	int total = 0;

	for (int i = 0; i < DIMENSIONS; i++) {
		total+= x[i] * x[i];
	}

	return total;
}

