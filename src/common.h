#ifndef COMMON_H_
#define COMMON_H_

enum {LOG_OBJECTIVE, LOG_RANDOM, LOG_STDOUT, LOG_SCALAR_ATRIBUTES, LOG_VECTOR_ATRIBUTES};

enum {ROSENBROOK, SPHERE, SCHWEFEL, ACKLEY, RASTRINGIN, GRIEWANK, SHUBER};

double sphere(double x[], int dimensions);
double rastringin (double solution[], int dimensions);
double griewank (double solution[], int dimensions);
double ackley (double solution[], int dimensions);
double rosenbrock (double solution[], int dimensions);
double schwefel(double solution[], int dimensions);
void my_seed(void);
double my_rand(int, int);
void logger(int destination, char *fmt, ...);
void allocate_resources();
void deallocate_resources();
double shuber (double solution[], int dimensions);


#endif
