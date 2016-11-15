#ifndef COMMON_H_
#define COMMON_H_

double sphere(double x[], int dimensions);
double rastringin (double solution[], int dimensions);
double griewank (double solution[], int dimensions);
double ackley (double solution[], int dimensions);
double rosenbrock (double solution[], int dimensions);
double schwefel(double solution[], int dimensions);
void my_seed(void);
double my_rand(int, int);
#endif
