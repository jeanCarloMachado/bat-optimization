bat metaheuristic
================

Implementations of the bat algorithm, with and without GPU.

Generating results
------------------

0 - ROSENBROOK
1 - SPHERE
2 - SCHWEFEL
3 - ACKLEY
4 - RASTRINGIN
5 - GRIEWANK
6 - SHUBER


Run CPU
-------

```sh
ITERATIONS=1000 BATS_COUNT=128 EVALUATION_FUNCTION=0 ./bat

```

Get averages
------------

```sh
RUN_TIMES=3 ITERATIONS=1000 BATS_COUNT=128 EVALUATION_FUNCTION=0 ./statistics GPU
```
