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

Quanto mais morcegos mais rápido o algorítmo converge para um
valor não ideal.


### O que acontece sem busca local?
O melhor converge para um valor longe do ideal e a média nunca
converge.

### O que acontece se fizer sempre a busca local?

Não existe grande diferença entre o melhor e o pior indivíduo e
o fitness fica próximo do número de dimensões.

### O que acontece se fazer a busca local de acordo com o
pulse-rate (forma padrão)?

A média fica mais longe do melhor mas o fitness fica próximo do
número de dimensões.



Valor centralizando em zero
---------------------------

Será que ele fica no centro pois o centro é a média de valores?
    Não, mudando os boundries para o domínio positivo os valores
    ainda ficam no centro.


Será o que problema da centralização da busca não está na busca local?
----------------------------------------------------------------------

Se a busca local é removida o algorítmo não chega nem perto de
convergir.



Discrepâncias com a implementação das proteinas
-----------------------------------------------

Processamentos repetidos dentro de loops.
A frequencia é gerada para cada uma das dimensões da posiçõ
temporária
Existe um pulse rate inicial para cada morcego.
Existem upper e lower boundries que não condizem com a literatura


BETA é (0,1) mas no código está (-1,1)

----

Seria interessante medir a proporção de posicões com 1 em relação
a proporção de posições com zero. Acho que o problema do algorítmo
possa ser que a busca local não é otimizada que chega.



Falhas do algorítmo
-------------------
Me parece certo dizer que o algorítmo proposto tem um falha de
fugir dos boundries com muita frequência, já aconteceu em 18% dos
números em uma execução.
