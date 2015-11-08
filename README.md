# WD_rotation

</Rotation of axially symetrical body in hydrostatic equillibrium. >

Unicode.cpp is the main numerical integrator. It contains three Runge-Kutta fourth order integrators by itself, but also it calculates Chandrasekhar limit, Stellar mass of White dwarfs, Radiuses (ekvatorial, polar and 45 degrees), density distribution over any direction in the stellar interior etc. 
It is used to obtain solution for the rotating equillibrium problem.

We assumed that the matter inside the white dwarf is in a state of relativistic completely degenerate electron gas and that the gravity is Newtonian. Solution to the equation of white dwarf structure has been approximated in the first order of perturbative expansion in perturbative parameter which depends on angular velocity and central density (Chandrasekhar-Milne expansion). 

Equation of state is already implemented inside the program. And it is higly relativistic completely degenerate gase. 
For wider use of this integrator, (note to self) it would be higly recommended to implement polytropic equation of state for its general properties.

Unicode.cpp by its self is designed as three main fields:

1) Constants
2) Integrator - The working machine 
3) Stellar properties

The working machine consists of three Fourth order Runge-Kutta formulas. They solve three differential equations. Solutions of those differential equations, obtained here, are for later use in the main part of the program.

Stellar properties calculates basic parameters of white dwarfs, using two "for"loops. One is for iterating central density of white dwarfs, and the outer loop shifts omega, the angular velocity of the star.

###

There is still much to be done in this code
