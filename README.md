### Code to accompany *[Quantum steering: a short review with focus on semi-definite programming](http://arxiv.org/abs/1604.00501)*
#### Daniel Cavalcanti and Paul Skrzypczyk

This repository provides a small collection of code which implements many of the semidefinite programs presented in the review article "*[Quantum steering: a short review with focus on semi-definite programming](http://arxiv.org/abs/1604.00501)*".

All code is written in MATLAB and requires:
- [CVX](http://cvxr.com/) - a Matlab-based convex modeling framework
- [QETLAB](http://www.qetlab.com/) - A MATLAB Toolbox for Quantum Entanglement

It has been tested on Matlab R2014a, and CVX 2.1 

The code comprises the following:

- Membership, helper, and misc:
  - [NSAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/NSAssemblage.m): determine whether a bipartite assemblage is a valid non-signalling assemblage or not<sup>§</sup>.
  - [LHSAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/LHSAssemblage.m): determine whether a biparitte assemblage has an LHS model or not<sup>§</sup>.
  - [genAssemblage](https://github.com/paulskrzypczyk/steeringreview/blob/master/genAssemblage.m): generate an assemblage starting from a quantum state and a set of measurements.
  - [validPOVMs](https://github.com/paulskrzypczyk/steeringreview/blob/master/validPOVMs.m): determine whether a set of POVMs is valid or not<sup>§</sup>.
  - [JMPOVMs](https://github.com/paulskrzypczyk/steeringreview/blob/master/JMPOVMs.m): determine whether a set of measurements is jointly measurable or not<sup>§</sup>.
  - [genRandProjMeas](https://github.com/paulskrzypczyk/steeringreview/blob/master/genRandProjMeas.m): generate a random set of projective measurements.
  - [genSinglePartyArray](https://github.com/paulskrzypczyk/steeringreview/blob/master/genSinglePartyArray.m): generate the single-party determinstic probability distributions
  - [bestSteeringMeasurements](https://github.com/paulskrzypczyk/steeringreview/blob/master/bestSteeringMeasurements.m): find the optimal measurements given a state and a steering functional.
  - [bestSteeringState](https://github.com/paulskrzypczyk/steeringreview/blob/master/bestSteeringState.m): find the optimal state given a set of measurements and a steering functional.
  - [findRadiusPolytopeInBlochSphere](https://github.com/paulskrzypczyk/steeringreview/blob/master/findRadiusPolytopeInBlochSphere.m): determine the radius of the largest ball which can fit inside a polytope contained in the Bloch ball determined by a set of measurements<sup>¶</sup>.
