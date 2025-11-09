# Structured finite volume (FVSTRCTRD) # {#nmFVSTRCTRD}

[TOC]

## Transformation to curvilinear coordinates
The transformation of the governing equations (NSE) to curvlinear coordinates

\f{align}{
  \tau = t, \qquad \xi = \xi(x,y,z,t), \qquad \eta = \eta(x,y,z,t), \qquad \zeta = \zeta(x,y,z,t)
\f}

yields the following transformed Navier-Stokes equations in conservative form

\f{align}{
\newcommand{\ddelt      } [1] {\frac{\partial #1}{\partial t}}
\newcommand{\ddelxi     } [1] {\frac{\partial #1}{\partial \xi}}
\newcommand{\ddeleta    } [1] {\frac{\partial #1}{\partial \eta}}
\newcommand{\ddelzeta   } [1] {\frac{\partial #1}{\partial \zeta}}
  \ddelt{\mv{\hat{Q}}} + \ddelxi{\mv{\hat{E}}} + \ddeleta{\mv{\hat{F}}} + \ddelzeta{\mv{\hat{G}}} = 0
\f}

where

\f{align}{
   \mv{\hat{Q}} = J\mv{Q}
\f}

and

\f{align}{
  \mv{\hat{E}} &= J\left(\xi_t\mv{Q} + \xi_x\mv{E} + \xi_y\mv{F} + \xi_z\mv{G}\right)\\
  \mv{\hat{F}} &= J\left(\eta_t\mv{Q} + \eta_x\mv{E} + \eta_y\mv{F} + \eta_z\mv{G}\right)\\
  \mv{\hat{G}} &= J\left(\zeta_t\mv{Q} + \zeta_x\mv{E} + \zeta_y\mv{F} + \zeta_z\mv{G}\right)
\f}

The quantities \f$\xi_i\f$, \f$\eta_i\f$\, and \f$\zeta_i\f$ denote metric terms and \f$J\f$ denotes the Jacobian.

To compute the transformed fluxes \f$\mv{\hat{E}}\f$, \f$\mv{\hat{F}}\f$, and \f$\mv{\hat{G}}\f$ with second-order accuracy central differences around the cell center are computed via
\f{align}{
    \left.\frac{\partial \mv{\hat{E}}}{\partial \xi}\right|_{i,j,k} &= \mv{\hat{E}}_{i+\frac{1}{2},j,k} - \mv{\hat{E}}_{i-\frac{1}{2},j,k}\\
    \left.\frac{\partial \mv{\hat{F}}}{\partial \eta}\right|_{i,j,k} &= \mv{\hat{F}}_{i,j+\frac{1}{2},k} - \mv{\hat{F}}_{i,j-\frac{1}{2},k}\\
\left.\frac{\partial \mv{\hat{G}}}{\partial \zeta}\right|_{i,j,k} &= \mv{\hat{G}}_{i,j,k+\frac{1}{2}} - \mv{\hat{G}}_{i,j,k-\frac{1}{2}}\,. 
\f}


## Convective flux discretization

The convective flux is computed with the Advection Upstream Splitting Method (AUSM) by [[LiouSteffen1993]]. In this method the inviscid flux is split into a convective part, computed by an upwind formulation, and a pressure part, computed by central differences. The inviscid flux at a cell interface, i.e., the surface that connects two cell volumes, by

\f{align}{
\newcommand{\vecfive   } [5] {\begin{pmatrix} #1\\ #2\\ #3\\ #4\\ #5 \end{pmatrix}}
  \mv{\hat{E}}_{inv} =
  J\Big[
  \vecfive{\rho U}{\rho Uu}{\rho Uv}{\rho Uw}{\rho U(E + p/\rho)} +
  \vecfive{0}{p\xi_x}{p\xi_y}{p\xi_z}{p\xi_t}\Big]
\f}

where the contravariant velocity is given by

\f{align}{
  U = u\xi_x + v\xi_y + w\xi_z + \xi_{t}\,.
\f}

The introduction of the speed of sound \f$a\f$ yields
\f{align}{
\newcommand{\vecfive   } [5] {\begin{pmatrix} #1\\ #2\\ #3\\ #4\\ #5 \end{pmatrix}}
  \mv{\hat{E}}_{L/R} =
\underbrace{
  \frac{U}{\left|\nabla \xi \right|a}
  }_{M_{L/R}}
  \underbrace{
  J\left|\nabla \xi\right|
  \vecfive{\rho a}{\rho au}{\rho av}{\rho aw}{\rho a\left(E + p/\rho\right)}
  }_{\mv{\hat{E}}^c_{L/R}}
  +
  \underbrace{
  J\vecfive{0}{p\xi_x}{p\xi_y}{p\xi_z}{p\xi_t}
  }_{\mv{\hat{E}}^p_{L/R}}
\f}

where \f$L\f$ and \f$R\f$ denote the left and right state of the surface. That is, the left or right state the interface Mach number

\f{align}{
  M_{L/R} = \frac{1}{2}\left[M_L + M_R\right]
\f}

where

\f{align}{
  M_{L,R} = \frac{U_{L,R}}{\left|\nabla \xi \right|a_{L,R}}
\f}

and

\f{align}{
  a_{L,R} = \sqrt{\gamma\frac{p_{L,R}}{\rho_{L,R}}}
\f}

Whether the left or right state is chosen for convection depends on the Mach number, i.e.,

\f{align}{
  \mv{\hat{E}}^c_{L/R} = \left\{
  \begin{array}{ll}
    (.)_L, &M_{L/R} \geq 0\\
    (.)_R, &M_{L/R} < 0    
  \end{array}
\right.\,.
\f}

As an example, the total discretized inviscid flux in the \f$\xi\f$-direction is computed by
\f{align}{
  \mv{\hat{E}}_{L/R} = \frac{1}{2}\left[M_{L/R}\left(\mv{\hat{E}}^c_L + \mv{\hat{E}}^c_R\right) +
  |M_{L/R}| \left(\mv{\hat{E}}^c_L - \mv{\hat{E}}^c_R\right)\right] + \mv{\hat{E}}^p_{L/R}
\f}

where

\f{align}{
  \mv{\hat{E}}^p_{L/R} = \frac{1}{2}\left[p_L + p_R\right]J\vecfive{0}{\xi_x}{\xi_y}{\xi_z}{\xi_t}\,.
  \f}

is the pressure part of the flux. To achieve a second-order accuracy the left  and right states are computed by a monotonic upstream-centered scheme for conservation laws (MUSCL) [[vanLeer1979]], such that

\f{align}{
\newcommand{\im      } [0] {{i - \frac{1}{2}}}
\newcommand{\ip      } [0] {{i + \frac{1}{2}}}
  \left.\mv{P}_L\right|_{i+\frac{1}{2}} &=
  \mv{P}_i + \frac{l_\ip}{\left(l_\ip + l_\im\right)^2}
  \left[l_\im \left(\mv{P}_{i+1} -\mv{P}_{i}\right) + l_\ip\left(\mv{P}_{i} -\mv{P}_{i-1} \right)\right]\\
  \left.\mv{P}_R\right|_{i+\frac{1}{2}} &=
  \mv{P}_{i-1} - \frac{l_\ip}{\left(l_\ip + l_{i+\frac{3}{2}}\right)^2}                                            
  \left[l_\ip \Delta_{i+\frac{3}{2}} \left(\mv{P}_{i+2} - \mv{P}_{i+1}\right) + l_{i+\frac{3}{2}}\left(\mv{P}_{i+1} - \mv{P}_{i}\right)\right]
  \f}

where

\f{align}{
l_{i\pm\frac{1}{2}} = \sqrt{\left(\pm\mv{x}_{i\pm 1} \mp
    \mv{x}_{i}\right)\cdot\left(\pm\mv{x}_{i\pm 1} \mp
    \mv{x}_{i}\right)}
\f}

is the distance function between the cell centers of neighbour cells in the chosen direction.

## Viscous flux discretization

The viscous flux is computed by a modified cell-vertex (MCV) scheme [Meinke1993], where the flux through each interface, e.g., \f$ \left(i+\frac{1}{2},j,k\right)\f$ and \f$ \left(i-\frac{1}{2},j,k\right)\f$ is determined by an average of the flux on the four cell interface corner points. For instance, for the \f$\mv{E}\f$-flux
\f{align}{
    \left.\frac{\partial\mv{\hat{E}}_{vis}}{\partial \xi}\right|_{i,j,k} =
    &\frac{1}{4}\left(
      \mv{\hat{E}}_{vis,i+\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2}} +
      \mv{\hat{E}}_{vis,i+\frac{1}{2},j-\frac{1}{2},k+\frac{1}{2}} +
      \mv{\hat{E}}_{vis,i+\frac{1}{2},j+\frac{1}{2},k-\frac{1}{2}} +
      \mv{\hat{E}}_{vis,i+\frac{1}{2},j-\frac{1}{2},k-\frac{1}{2}}
    \right)\\
    - &\frac{1}{4}\left(
      \mv{\hat{E}}_{vis,i-\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2}} +
      \mv{\hat{E}}_{vis,i-\frac{1}{2},j-\frac{1}{2},k+\frac{1}{2}} +
      \mv{\hat{E}}_{vis,i-\frac{1}{2},j+\frac{1}{2},k-\frac{1}{2}} +
      \mv{\hat{E}}_{vis,i-\frac{1}{2},j-\frac{1}{2},k-\frac{1}{2}}
    \right)
\f}

The corner fluxes are computed by central differences using the variables on the surrounding eight cell centers, e.g., for \f$\frac{\partial u}{\partial \xi}\f$ at the corner point \f$ \left(i+\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2}\right)\f$

\f{align}{
    \left.u_{\xi}\right|_{i+\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2}} =&\frac{1}{4}\left(u_{i+1,j,k} + u_{i+1,j+1,k} + u_{i+1,j,k+1} + u_{i+1,j+1,k+1} \right) -\\
    &\frac{1}{4}\left(u_{i,j,k} + u_{i,j+1,k} + u_{i,j,k+1} + u_{i,j+1,k+1}\right)\,.
\f}

## Temporal integration

Combining the discretized spatial derivatives into the discrete right-hand-side operator \f$RHS\left(t;\cdot \right)\f$ the semi-discrete equation
\f{align}{
  \frac{\partial \mv{\hat{Q}}}{\partial t} = RHS\left(t;\mv{\hat{Q}}\right)
\f}

is then integrated in time by a Runge-Kutta scheme, such that

\f{align}{
  \begin{array}{ccl}
    \mv{\hat{Q}}^{(0)} &=& \mv{\hat{Q}}^n\\
    &\vdots& \\
    \mv{\hat{Q}}^{(k)} &=& \mv{\hat{Q}}^{(0)} + \alpha_k\Delta t  RHS\left(t^n + \alpha_{k-1}\Delta t; \mv{\hat{Q}}^{(k-1)}\right)\quad k = 1,\dots ,s\,,\\
                      &\vdots& \\
    \mv{\hat{Q}}^{n+1} &=& \mv{\hat{Q}}^{(k)}
    
  \end{array}
  \f}

In most computations the coefficients are set to \f$\alpha_k = (1/4, 1/6, 3/8, 1/2, 1)\f$ such that a second-order accurarcy in time is achieved. 

## References
* M.-S. Liou, C. Steffen, A new flux splitting scheme, J. Comput. Phys. 107 (1993)
23–39. [10.1016/j.compfluid.2014.07.004][LiouSteffen1993]
* B. van Leer, Towards the ultimate conservative difference scheme. V. A second-
order sequel to Godunov’s method, J. Comput. Phys. 32 (1979) 101–136. [10.1016/0021-9991(79)90145-1][vanLeer1979]
* M. Meinke, Numerische Lösung der Navier-Stokes-Gleichungen für instationäre Strö-
mungen mit Hilfe der Mehrgittermethode, Ph.D. thesis, Institute of Aerodynamics,
RWTH Aachen University, 1993

[LiouSteffen1993]: https://doi.org/10.1016/j.compfluid.2014.07.004
[vanLeer1979]: https://doi.org/10.1016/0021-9991(79)90145-1