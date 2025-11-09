# Lattice Boltzmann equation (LBE) # {#mmLBE}
In the kinetic gas theory molecules and atoms are modelled with hard, point-like spheres,
that collide ellastically with other molecules or the domain boundary.
Based on the prodcedure in [Hänel],
the mean molecular velocity for one of these spheres \f$ \xi_0=\sqrt{3RT} \f$ can be derived,
where \f$ R \f$ is the specific gas constant and \f$ T \f$ the temperature. 
To draw conclusions about the macroscopic behaviour of the gas 
one would have to calculate the trajectory of each particle,
which is not possible in most cases, due to the large numbe rof particles.
Therefore, a statistical description is used to determine the macroscopic variables of the gas.

## Particle probability distribution function (PPDF)
The particle probability distribution function (PPDF) \f$ f(\mv{x},\mv{\xi},t) \f$
represents the probability of particles to appear in momentum space at the
position \f$ \mv{x} \f$ and time \f$ t \f$ with a velocity of \f$ \mv{\xi} \f$.

Macroscopic quantities such as density, velocity, or energy
(\f$\rho , \mv{u}, E\f$) are connected to this mesoscopic quantity \f$ f \f$ by
its raw moments
\f{align}{
  \rho(\mv{x},t)            &= \iiint               f(\mv{x},\mv{\xi},t) d\mv{\xi}  \,, \\
  \rho(\mv{x},t)\mv{u}(\mv{x},t) &= \iiint  \mv{\xi}     f(\mv{x},\mv{\xi},t) d\mv{\xi}  \,, \text{and}\\
  2\rho(\mv{x},t)E(\mv{x},t) &= \iiint |\mv{\xi}|^2  f(\mv{x},\mv{\xi},t) d\mv{\xi} \,.
\f}

However, to calculate the moments of the PPDFs and thus the macroscopic flow variables,
the PPDFs must be determined first.
The transport equation of the PPDF is the Boltzmann equation.

## Boltzmann equation
The Boltzmann equation describes the temporal evolution of this PPDF reading
\f{equation}{
  \fracpart{f}{t} + \mv{\xi} \cdot \nabla_x{f} + \frac{\mv{F}}{m} \nabla_\xi{f} = \Omega(f) \,,
\f}
with the force vector \f$ F \f$, the molecule mass \f$ m \f$,
and the collision operator \f$ \Omega \f$ as its source term accounting for the
effect of the momentum exchange of particles colliding with each other.
Although solutions of the Boltzmann equation exist for special conditions,
it is usually not possible to solve the Boltzmann equation under general conditions.
Hence, Bhatnagar, Gross, and Krook (BGK) developed a simplified model of the complex collision term [Bhatnagar1954]
\f{equation}{
  \fracpart{f}{t} + \mv{\xi} \cdot \nabla_x{f} + \frac{\mv{F}}{m} \nabla_\xi{f} = \omega(f^{eq}-f) \,,
\f}
with the Maxwell distribution function
\f{equation}{
  f^{eq}(\mv{x}, \mv{\xi}, t)=\frac{\rho(\mv{x},t)}{(2 \pi RT(\mv{x},t))^{\frac{2}{3}}} \cdot e^{-\frac{(\mv{\xi}-\mv{u}(\mv{x}, t))^2}{2RT(\mv{x},t)}} \,.
\f}

##Forcing

TODO: Forcing

## References

* D. Hänel, Molekulare Gasdynamik, Einführung in die kinetische Theorie der
Gase und Lattice-Boltzmann-Methoden.

* P. L. Bhatnagar, E. P. Gross, and M. Krook, A Model for Collision Processes
in Gases. I. Small Amplitude Processes in Charged and Neutral One-Component
Systems, Phys. Rev., 94, 511–525, May 1954, [10.1103/PhysRev.94.511][Bhatnagar1954].

[Bhatnagar1954]: https://doi.org/10.1103/PhysRev.94.511

