bat-optimization
================

Implementations of the bat algorithm, with and whitout GPU.


Rationing
---------


O schewel tende a se concentrar no boundry mas converge.

O Rosenbrook não converge, e a posição tende a ser o zero.
O valor do rosenbrook tente a ficar próximo ao número de
dimensões.

### Os números tendem mais a zero?

Fiz um contador e os números tendiam mais a 1 no rosenbrook. Porém
eles eram perdidos até o final da execução. Talvez a probabilidade
de aceitar números ruins esteja muito alta no final da execução.


Hipóteses 
----------

Será que o algorítmo rosenbrook não converge pois não funciona
para o problema? 



### Sat Nov 19 16:31:27 BRST 2016

The values are not converging anymore I'll have to see why this is
happening.

The pulse rate seems to decay the same way of the convergence,
maybe if I make it endure longer the convergence will follow.


### Sun Nov 20 11:13:01 BRST 2016

It seems that the algorithms always converges on the 20th
iteration I should discover why this happens.
