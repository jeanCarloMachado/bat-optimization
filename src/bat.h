#ifndef BAT_H_
#define BAT_H_


enum {LOG_OBJECTIVE, LOG_RANDOM, LOG_STDOUT, LOG_SCALAR_ATRIBUTES, LOG_VECTOR_ATRIBUTES};

enum {ROSENBROOK, SPHERE, SCHWEFEL, ACKLEY, RASTRINGIN, GRIEWANK, SHUBER};

void logger(int destination, char *fmt, ...);
char* get_function_name(int index);

#define DUMP_DIR "./dump"

#endif
