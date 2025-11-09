# Ffowcs Williams and Hawkings (FWH) equation # {#mmFWH}
The Ffowcs Williams and Hawkings (FWH) equation [[WilliamHawkings1969]] is
obtained by an exact rearrangement of the [Navier-Stokes equations](@ref mmNS),
i.e, the conservation of mass, momentum, and energy, into an inhomogeneous wave
equation featuring a non-linear right-hand side
\f{align}{
  \left( \fracpart{^2}{t^2} - c_s^2 \fracpart{^2}{x_i^2} \right) (H(s) \rhop)
  &=
  \fracpart {^2}{x_i x_j}(T_{ij} H(s)) 
  - \fracpart{}{x_i} (F_i \delta(s)) 
  + \fracpart{}{t} (Q \delta(s)) &&\\
  T_{ij}  &= \rho u_i u_j + P_{ij} - c_s^2 \rhop \delta_{ij}          &&\text{(Quadrupole term)}  \\
  F_i     &= (P_{ij} + \rho u_i (u_j - v_j)) \fracpart{s}{x_j}        &&\text{(Dipole term)}      \\
  Q       &= ( \rho_\infty v_i + \rho(u_i - v_i) ) \fracpart{s}{x_i}  &&\text{(Monopole term)} \,,
\f}
where \f$\delta_{ij}\f$ is the Kronecker delta, \f$\mm{P}=p\delta_{ij}\f$ is the
compressive stress tensor, \f$\rhop\f$ is the perturbed density, and \f$H(\cdot)\f$ is 
the Heaviside step function.
The function \f$s\f$ is a support function, which takes zero value on a surface
surrounding the volume of acoustic source terms \f$V_a\f$, negative value in its
inside, and positive values outside of this volume. From now on this surface is
named FWH permeable surface.
With \f$H()\f$ being the Heaviside step function it follows that
\f{equation*}{
  H(s(\mv{x}, t)) = 
  \begin{cases}
    1,  & \mv{x} \notin V_a \\
    0,  & \mv{x} \in    V_a
  \end{cases}
  \,.
\f}
The FWH permeable surface \f$ s=0 \f$ is moving with the velocity \f$ \mv{v}
\f$.

## Transformation into frequency space
Following [[Lockard2000]], the FWH equation can be solved more efficiently in
frequency domain than in time domain. Therefore, a uniform rectilinear motion of
the surface \f$ s=0 \f$ with the speed \f$ \mv{U} \f$ needs to be assumed.
Applying a Galilean transformation
\f{equation*}{
  \mv{y}            = \mv{x} + \mv{U} t  , \qquad
  \fracpart{}{x_i}  = \fracpart{}{y_i}   , \qquad
  \fracpart{}{t}    = \fracpart{}{t'} + U_i \fracpart{}{y_i}
\f}
and a Fourier transformation, the final governing equation consists of two
surface integrals and one volume integral reading
\f{aligned}{
  H(s)c_s^2 \rhop(\mv{y_o}, \omega) 
  = 
  &- \oint_{s=0} \iu \omega Q(\mv{y_s}, \omega) G(\mv{y_o}; \mv{y_s}) dA \\
  &- \oint_{s=0} F_i(\mv{y_s}, \omega) \fracpart{G(\mv{y_o}; \mv{y_s})}{y_{s,i}} dA \\
  &- \int_{s>0} T_{ij}(\mv{y_s}, \omega) H(s) \fracpart{^2G(\mv{y_o}; \mv{y_s})}{y_{s,i} \partial y_{s,j}} \mv{dy_s}
  \,,
\f}
with the angular frequency \f$\omega\f$ and the Green's function
\f$G(\cdot;\cdot)\f$. Source and observer coordinates are given by
\f$\mv{y_s}\f$ and \f$\mv{y_o}\f$, respectively. Here, observer refers to
location of interests where the acoustic signal is to be investigated and source
describes all locations contributing to the resulting signal.

As the term \f$ H(s(\mv{x}))=0 , \forall \mv{x} \in V_A \f$ the sum of the three
integrals must zero inside of the volume of acoustic source terms per
definition.

## Simplification
In it's integral form the FWH equation features two surface integrals which are
to be solved on the chosen FWH permeable surface and one volume integral to be
solved only in the volume that is not enclosed by the FWH permeable surface.
Thus, if all significant contributions of acoustic quadrupole sources
represented by the Lighthill stress tensor \f$ T_{ij} \f$ being part of the
volume integral are contained within the enclosed volume it is a good and common
approximation to neglect the calculation of the volume integral.

## References
* Williams, J. E. F., and Hawkings, D. L., “Sound generation by turbulence and surfaces in arbitrary motion,” Phil. Trans. R. Soc., London A, 1969, pp. 321–342. [jstor.org/stable/73790][WilliamHawkings1969]
* Lockard, D. P., “An efficient, two-dimensional implementation of the Ffowcs Williams and Hawkings equation,” J. Sound Vibr.,
Vol. 229, No. 4, 2000, pp. 897–911. [doi.org/10.1006/jsvi.1999.2522][Lockard2000].

[WilliamHawkings1969]: https://www.jstor.org/stable/73790
[Lockard2000]: https://doi.org/10.1006/jsvi.1999.2522
