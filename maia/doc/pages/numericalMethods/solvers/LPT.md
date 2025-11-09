# Lagrangian particle tracking (LPT) # {#nmLPT}

A modified second-order accurate Euler method is used to solve 
the particle motion equation (\ref{partveq}). 
The differential equations for particle temperature- and 
mass change (Eqs.~\ref{parteqm} and~\ref{parteqt}) are 
solved by an iterative implicit Euler method.
Continuous phase state variables (\f$\rho, \mv{u}, T, p, Y\f$) are 
interpolated first-order accurately to the parcel position. 
A stencil containing all surrounding cells of the parcel and a
distance weight is applied.
The same weights are used for the redistribution of the parcel source 
terms for the continuous phase description.
Particle-wall collisions are modeled as hard-sphere collisions based on 
boundary surface normal and position.
