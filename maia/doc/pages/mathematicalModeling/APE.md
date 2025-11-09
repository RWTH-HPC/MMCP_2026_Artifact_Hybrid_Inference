# Acoustic pertubation equation (APE) # {#mmAPE}

Source: [Schlottke-Lakemper2019][Schlottke-Lakemper2019]

The acoustic perturbation equations (APE) are derived from the conservation equations and
modified to retain only acoustic modes (@ref EwertSchroeder03). The APE-4 system
can be written in non-dimensional form as
\f{align}{
  \frac{\partial \mv{u}'}{\partial t} + \mv{\nabla} \left( \bar{\mv{u}} \cdot \mv{u}' +
  \frac{p'}{\bar{\rho}} \right) &= \mv{q}_m, \\
  \frac{\partial p'}{\partial t} + \bar{c}^2 \mv{\nabla} \cdot \left( \bar{\rho}\mv{u}' +
  \bar{\mv{u}} \frac{p'}{\bar{c}^2}\right) &= 0.
\f}
The unknowns of the APE are acoustic velocity and pressure fluctuations, which are denoted by prime \f$(\cdot)'\f$. They are
defined by \f$\phi' := \phi - \bar{\phi}\f$, where the bar \f$\overline{(\cdot)}\f$ indicates mean
quantities that have to be determined prior to the CAA simulation by time-averaging the flow
solution.
In the momentum equations, only the contribution of the Lamb vector source is usually considered for problems with noise generation dominated by vorticity, which
is given by
\begin{equation}
  \mv{q}_m = - (\mv{\omega} \times \mv{u})',\label{eqn:ape_q_m}\\
\end{equation}
with the vorticity vector \f$\mv{\omega} = \nabla \times \mv{u}\f$.
It is calculated using data from the flow simulation, i.e., the perturbed Lamb vector
in \f$\eqref{eqn:ape_q_m}\f$ is based on hydrodynamic quantities.
To discretize the APE with a discontinuous Galerkin scheme, they are rewritten in conservative form.
Based on the formulation in [Schlottke-Lakemper2017][Schlottke-Lakemper2017], the system reads

\f{align}{
  \frac{\partial \mv{u}'}{\partial t} + \mv{\nabla} \left( \bar{\mv{u}} \cdot \mv{u}' +
  \frac{p'}{\bar{\rho}} \right) &= \mv{q}_m, \\
  \frac{\partial p'}{\partial t} + \mv{\nabla} \cdot \left( \bar{c}^2\bar{\rho}\mv{u}' +
    \bar{\mv{u}} p'\right) - \underbrace{\left(\bar{\rho}\mv{u}' + \bar{\mv{u}}
  \frac{p'}{\bar{c}^2} \right) \cdot \nabla \bar{c}^2}_{s_\text{cons}} &= 0,
\f}
where \f$s_\text{cons}\f$ is an additional term introduced by reformulating the APE-4 system.


### References:

* M. Schlottke-Lakemper, A. Niemoeller, M. Meinke, W. Schroeder,
Efficient parallelization for volume-coupled multiphysics simulations on hierarchical Cartesian grids, Computer Methods in Applied Mechanics and Engineering, Volume 352, pp 461-487, 2019, [10.1016/j.cma.2019.04.032][Schlottke-Lakemper2019].
[Schlottke-Lakemper2019]: https://doi.org/10.1016/j.cma.2019.04.032

* R. Ewert, W. Schroeder, Acoustic perturbation equations based on flow decomposition via source filtering, J. Comput. Phys. 188 (2003) 365–398, [10.1016/S0021-9991(03)00168-2][EwertSchroeder03].
[EwertSchroeder03]: http://dx.doi.org/10.1016/S0021-9991(03)00168-2

* M. Schlottke-Lakemper, H. Yu, Sven Berger, M. Meinke, W. Schroeder, A fully coupled hybrid computational aeroacoustics method
on hierarchical Cartesian meshes, Comput. Fluids 144 (2017) 137–153, [10.1016/j.compfluid.2016.12.001][Schlottke-Lakemper2017].
[Schlottke-Lakemper2017]: http://dx.doi.org/10.1016/j.compfluid.2016.12.001

