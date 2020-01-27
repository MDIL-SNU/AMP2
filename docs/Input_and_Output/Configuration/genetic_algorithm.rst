genetic_algorithm
-----------------
AMP2 uses genetic algorithm to identify thr most stable magnetic spin ordering. The population in a generation is 90
and the number of generations is ten times of the number of magnetic atoms. In each step, the parents are selected
according to Ising energy of each gene. The fitness function (:math:`f`) for selection is

.. math::
    f = E_{\textrm{max}} - E + (E_{\textrm{max}} - E_{\textrm{min}}) \times \gamma / (\gamma - 1)

where :math:`\gamma` is election coefficient (2). Using small :math:`\gamma`, the parents are more randomly chosen.

Each population consists of 15 the most stable candidates, 30 the candidates made by crossover scheme, 
30 the candidates made by mutation scheme and 15 randomly generated candidates.

