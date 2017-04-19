
enum {LOG_OBJECTIVE, LOG_RANDOM, LOG_STDOUT, LOG_SCALAR_ATRIBUTES, LOG_VECTOR_ATRIBUTES};
enum {ROSENBROOK, SPHERE, SCHWEFEL, ACKLEY, RASTRINGIN, GRIEWANK, SHUBER};
double schwefel(double solution[], int dimensions);
void initialize_function(int evaluation_function);

