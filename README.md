A Monte Carlo algorithm written in C/C++, to calculate the effective
elastic coefficients of a membrane in thermal equilibrium with
its environment. The mesoscopic description by Helfrich and Canham is
employed.


The following cases are studied:

a) Untethered-free membrane.

b) Tethered membrane. Individual degrees of freedom are tethered after
   a random selection drawn out of a uniform probability distribution.

c) Block-tethered membrane. Regions/blocks of user-defined radius are tethered.
   In this case, interconnected tethered regions may arise.