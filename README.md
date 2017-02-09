bat-optimization
================

Implementations of the bat algorithm, with and whitout GPU.

numero aleatorio
pseudo codigo GPU
usar ackley, griewank, rastringin, rosenbrook
200 dim
(procurar limites)

Padrão de escrita: i3e max 10 pg

Modelo de carta de referência.
Usar boundries menores.

Event to submit
---------------
21 dezembro IEEE parallel distributing processsing PDCO 2017


Todo:

Figura da arquitetura - como as threads
Procurar tradeoffs de números aleatórios em GPU.

Colocar desvio padrão.

Gráfico de convergência/ CPU,GPU

Mandar o paper pro parpinelli.

Testar mil dimensões.

Validar a possibilidade de colocar mil dimensões.

SPAA - simposio da ACM

Generating results
------------------

0 - ROSENBROOK
1 - SPHERE
2 - SCHWEFEL
3 - ACKLEY
4 - RASTRINGIN
5 - GRIEWANK
6 - SHUBER


```
TOTAL=3 ITERATIONS=1000 BATS_COUNT=128 EVALUATION_FUNCTION=0 ./statistics GPU

```
