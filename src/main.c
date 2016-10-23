#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#define DIMENSIONS 2
#define MAX_ITERATIONS 700000
#define BATS_COUNT 40
#define FREQUENCY_MIN 0
#define FREQUENCY_MAX 100
#define LOUDNESS_MIN 1
#define LOUDNESS_MAX 100
#define ALFA 0.5
#define LAMBDA 0.1
#define DUMP_DIR "/home/jean/projects/bat-optimization/dump"

#define DEBUG 1
#define DEBUG_RANDOM 1

#define LOG_FILE_MAIN 1
#define LOG_FILE_RANDOM 2

int RUN_TIME;
FILE *LOG;
FILE *LOG_RANDOM;

struct bat {
	double pulse_rate; //or frequency
	double loudness;
	double frequency;
	double position[DIMENSIONS];
	double velocity[DIMENSIONS];
};


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
struct bat get_worst(struct bat bats[]);
struct bat get_average(struct bat bats[]);
void my_print(int destination, char *fmt, ...);
void allocate_resources(void);
void deallocate_resources();

int main() {
	allocate_resources();
	struct bat bats[BATS_COUNT];
	struct bat best;
	struct bat worst;
	struct bat average;
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

			if (my_random(0,1) < bats[j].loudness || objective_function(candidate) < objective_function(*current)) {
				memcpy(current->position, candidate.position, sizeof candidate.position);
				current->pulse_rate = 1 - exp(-LAMBDA*iteration);
				current->loudness =  ALFA*current->loudness;
			}

			best = get_best(bats);	

		}
	}

	my_print(LOG_FILE_MAIN, "BEST");
	best = get_best(bats);
	print_bat(best);
	my_print(LOG_FILE_MAIN, "AVERAGE");
	average = get_average(bats);
	print_bat(average);
	my_print(LOG_FILE_MAIN, "WORST");
	worst = get_worst(bats);
	print_bat(worst);


	deallocate_resources();
	return 0;
}

void allocate_resources()
{
	RUN_TIME = time(NULL);

	char fileName[100];
	sprintf(fileName, "%s/%i-main", DUMP_DIR, RUN_TIME);
	LOG = fopen(fileName,"w");
	if (LOG == NULL)
	{
		printf("Error opening file %s !\n", fileName);
		exit(1);
	}

	sprintf(fileName, "%s/%i-random", DUMP_DIR, RUN_TIME);
	LOG_RANDOM = fopen(fileName,"w");
	if (LOG_RANDOM == NULL)
	{
		printf("Error opening file %s !\n", fileName);
		exit(1);
	}



}

void deallocate_resources()
{
	fclose(LOG);
	fclose(LOG_RANDOM);
}

void my_print(int destination, char *fmt, ...)
{
	char formatted_string[6666];

	va_list argptr;
	va_start(argptr,fmt);
	vsprintf(formatted_string, fmt, argptr);
	va_end(argptr);

	/* printf("%s",formatted_string); */

	if (destination == LOG_FILE_MAIN) 
		fprintf(LOG,"%s",formatted_string);
	else if (destination == LOG_FILE_RANDOM)
		fprintf(LOG_RANDOM,"%s",formatted_string);
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
	my_print(LOG_FILE_MAIN, " = BAT = \n");
	my_print(LOG_FILE_MAIN, "\tFrequency: %f\n", bat.frequency);
	my_print(LOG_FILE_MAIN, "\tLoudness: %f\n", bat.loudness);
	my_print(LOG_FILE_MAIN, "\tPulse-rate: %f\n", bat.pulse_rate);

	my_print(LOG_FILE_MAIN, "\tVelocity:\n");
	for (int i = 0; i < DIMENSIONS; i++) {
		my_print(LOG_FILE_MAIN, "\t[%d] %f \n", i, bat.velocity[i]);
	}

	my_print(LOG_FILE_MAIN, "\tPosition:\n");
	for (int i = 0; i < DIMENSIONS; i++) {
		my_print(LOG_FILE_MAIN, "\t[%d] %f \n", i, bat.position[i]);
	}
}
struct bat get_worst(struct bat bats[])
{
	double current_worst_val; 
	double current_val;

	current_val = current_worst_val = objective_function(bats[0]);
	struct bat current_worst_bat = bats[0];
	for (int i = 0; i < BATS_COUNT; i++) {
		current_val = objective_function(bats[i]);
		if (current_val > current_worst_val) {
			current_worst_bat = bats[i];
			current_worst_val = current_val;
		}
	}

	return current_worst_bat;
}


struct bat get_average(struct bat bats[])
{
	struct bat average;

	for (int i = 0; i < BATS_COUNT; i++) {
		average.frequency =  bats[i].frequency;
		average.loudness =  bats[i].loudness;
		average.pulse_rate =  bats[i].pulse_rate;
		for (int j = 0; j < DIMENSIONS; j++) {
			average.position[j]+= bats[i].position[j];
			average.velocity[j]+= bats[i].velocity[j];
		}
	}


	average.frequency =  average.frequency / BATS_COUNT;
	average.loudness =  average.loudness / BATS_COUNT;
	average.pulse_rate =  average.pulse_rate / BATS_COUNT;
	for (int j = 0; j < DIMENSIONS; j++) {
		average.position[j] = average.position[j] / BATS_COUNT;
		average.velocity[j] = average.velocity[j] / BATS_COUNT;
	}

	return average;
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
	double result;
	if (inferior == 0 && superior == 1) {
		result = (double)rand() / (double)RAND_MAX ;

		my_print(LOG_FILE_RANDOM, "0-1: %f\n", result, "debug");
		return result;
	}

	result = rand () % (int) superior;
	my_print(LOG_FILE_RANDOM, "0-100: %f\n", result, "debug");
	return result;
}

void initialize_bats(struct bat bats[])
{
	for (int i = 0; i < BATS_COUNT; i ++ ) {
		bats[i].pulse_rate = 0;
		bats[i].frequency = 0;
		bats[i].loudness = 1;

		for (int j = 0; j < DIMENSIONS; j++) {
			bats[i].velocity[j] = 0;
			bats[i].position[j] = my_random(0, RAND_MAX);
		}
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

