### Code to accompany *[All entangled states can demonstrate non-classical teleportation](https://arxiv.org/abs/1607.03249)*  
#### Daniel Cavalcanti, Paul Skrzypczyk and Ivan Šupić

This repository provides a (very) small collection of code which implements the semidefinite program presented in the article  
*[All Entangled States can Demonstrate Nonclassical Teleportation](https://doi.org/10.1103/PhysRevLett.119.110501)*.  
Daniel Cavalcanti, Paul Skrzypczyk, and Ivan Šupić  
Phys. Rev. Lett. **119**, 110501 (2017)
>>>>>>> f91aebfb9a965131673e0c6d5cb84a8f98d37711

All code is written in MATLAB and requires:
- [CVX](http://cvxr.com/) - a Matlab-based convex modeling framework
- [QETLAB](http://www.qetlab.com/) - A MATLAB Toolbox for Quantum Entanglement

It has been tested on Matlab R2016b, and CVX 2.1 

The code comprises the following:

  - [genTeleportationData](https://github.com/paulskrzypczyk/nonclassicalteleportation/blob/master/genTeleportationData.m): generates the unnormalised states prepared for Bob, when Alice performs a joint measurement on the (unknown) input state and her half of a shared state.
  - [teleportationRandomRobustness](https://github.com/paulskrzypczyk/nonclassicalteleportation/blob/master/teleportationRandomRobustness.m): calculates the random teleportation robustness of a given set of teleportation data. 
  
Slides from a presentation about this work can be found [here](https://github.com/paulskrzypczyk/Talks/blob/master/Teleportation%20-%20TyQI%202017.pdf).
