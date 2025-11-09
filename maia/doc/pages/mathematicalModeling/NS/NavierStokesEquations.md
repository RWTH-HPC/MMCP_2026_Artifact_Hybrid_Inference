# Navier-Stokes Equations # {#mmNavierStokesEquations}

The flow of a compressible fluid is governed by the general
conservation equations for mass, momentum and energy in integral form.
In the following, all variables are assumed to be made nondimensional,
see section [Nondimensionalization of the Navier-Stokes Equations](@ref mmNonDim). The integral form reads

\begin{equation}
\int_\tau \fracpart{\mv{Q}}{t} d\tau + \oint_A \mv{H} \cdot \mv{n} dA =
\int_\tau \mv{F_{\tau}} \,\, d\tau 
\end{equation}

for an infinitesimal control volume \f$ \tau \f$, its surface \f$ A
\f$ and the normal vector \f$ \mv{n} \f$.  Herein, \f$
\mv{Q}=\left(\rho, \rho \mv{u}, \rho E\right)^T \f$ is the vector of
conservative variables defined by the density \f$ \rho \f$, the
velocity \f$ \mv{u} \f$, the total specific energy \f$
E=e+|\mv{u}|^2/2 \f$, and the specific internal energy \f$ e \f$, and
\f$ \mv{H} \f$ is the flux vector. \f$ \mv{F_{\tau}} \f$ denotes a
source term including, e.g., external forces exerted on the volume \f$ \tau
\f$ such as gravity.

The flux vector can be decomposed into an inviscid part \f$ \mv{H}^i
\f$ and a viscous part \f$ \mv{H}^v \f$:

\f{equation}
\mv{H} =
\mv{H}^i - \mv{H}^v =
\left( \begin{array}{c}
\rho\mv{u} \\
\rho\mv{u}\mv{u} + p \\
\mv{u}(\rho E + p)
\end{array} \right) +
\frac{1}{Re_0}
\left( \begin{array}{c}
0 \\
\mm{\tau} \\
\mm{\tau} \mv{u} + \mv{q}
\end{array}\right)
\, ,
\f}

where \f$ p \f$ denotes the pressure, \f$ \mm{\tau} \f$ the stress
tensor, and \f$ \mv{q} \f$ the vector of heat conduction.  The
Reynolds number \f$ Re_0=\rho^*_0 a^*_0 L^* /\mu^*_0 \f$ is defined using the
stagnation point values for the density \f$ \rho^*_0 \f$, the dynamic
viscosity \f$ \mu^*_0 \f$, the freestream reference velocity \f$ a^*_0
\f$, and the reference length \f$ L^* \f$, where the superscript \f$ ^* \f$ denotes  dimensional quantities.

In case of a Newtonian fluid with zero bulk viscosity, the stress
tensor can be written as

\begin{equation}
  \mm{\tau} = - 2\mu\mm{S} + \frac{2}{3} \mu (\mv{\nabla}\boldsymbol{\cdot}\mv{u}) \mm{I}
  \, ,
\end{equation}

where \f$ \mm{S} =
\left[\mv{\nabla}\mv{u}+(\mv{\nabla}\mv{u})^T\right]/2 \f$, i.e., the
strain rate tensor, and \f$ \mm{I} \f$ is the identity matrix. The
temperature dependence of the dynamic viscosity \f$ \mu \f$ can be
evaluated according to Sutherland's law

\begin{equation}
\frac{\mu^*}{\mu^*_{ref}} =
\left(\frac{T^*}{T^*_{ref}}\right)^{3/2} \frac{T^*_{ref} + T^*_S}{T^* + T^*_S}
\, ,
\end{equation}

introducing the reference viscosity \f$ \mu^*_{ref} \f$, the static
temperature \f$ T^* \f$, the Sutherland reference temperature \f$
T^*_{ref} \f$, and the Sutherland constant \f$ T^*_S \f$. For air,
these constants are usually chosen to be \f$ \mu^*_{ref} = 1.716
10^{−5} \frac{kg}{ms} \f$, \f$ T^*_{ref} = 273.15 \mathrm{K} \f$, \f$
T^*_S=110.4 \mathrm{K} \f$. In nondimensional form, Sutherland's law is
given by

\begin{equation}
\frac{\mu}{\mu_{ref}} = \left(\frac{T}{T_{ref}}\right)^{3/2} \frac{T_{ref} + T_S}{T + T_S}
\, .
\end{equation}

The Sutherland constant is defined within @maia as [m_sutherlandConstant](@ref m_sutherlandConstant) and used in the macro SUTHERLANDLAW.


The thermal conductivity is expressed by Fourier's law

\begin{equation}
\mv{q} = -\frac{k}{Pr(\gamma-1)} \mv{\nabla} T
\, .
\end{equation}

The Prandtl number \f$ Pr=\mu^*_0 c^*_p/k^*_0 \f$ contains the
specific heat at constant pressure \f$ c^*_p \f$, \f$ k^* \f$ is the
thermal conductivity, \f$ \gamma=c^*_p/c^*_v \f$ is the ratio of specific
heats, where \f$ c^*_v \f$ is the specific heat at constant volume.

Assuming a constant Prandtl number and constant specific heats, the temperature dependence of the thermal conductivity
becomes \f$ k(T)=\mu(T) \f$. Closure of the equations is achieved via the
caloric equation of state \f$ e^*=c^*_v T^* \f$ and the ideal gas law
\f$ p^*=\rho^* R^* T^* \f$, with \f$ R^*=c^*_p-c^*_v \f$ as the specific gas constant.

In nondimensional form the ideal gas law reads
\begin{equation}
T = \gamma \frac{p}{\rho} 
\, .
\end{equation}



# Reynolds Averaged Navier-Stokes Equations # {#mmRANS}

The Reynolds-averaged Navier-Stokes (RANS) equations are obtained after time averaging the Navier-Stokes equations over all turbulent scales, where the time average for any quantity \f$ \phi \f$ is computed as

\begin{equation}
\mta{\phi} = \int_{\Delta t} \phi (t) \, dt
\quad ,
\end{equation}	
where \f$ \Delta t \f$ denotes a time intervall which is large enough to filter out all turbulent scales.

With the Reynolds decomposition the time averaged solution of the flow field
variables such as the time averaged velocity \f$ \mta{u} \f$, can be decoupled
from the time-varying, fluctuating component \f$ \up \f$, which accounts
for turbulent fluctuations in the flow field. In cases where
compressibility of the fluid has to be considered, a Favre
(density-weighted) decomposition of the Navier-Stokes equations is conducted to avoid additional density-velocity correlations. The
Favre average is denoted by a \f$ \mfa{\,} \f$ and is defined as \f$
\mfa{u} = \mta{\rho u} / \mta{\rho} \f$. Favre averaging is applied
to all variables except the density \f$ \rho \f$ and pressure \f$ p \f$.

The Reynolds decomposition is defined as
\begin{eqnarray}
  \mv{u}(\mv{x},t) = \mta{\mv{u}}(\mv{x}) + \mv{\up}(\mv{x},t)
  \quad ,
\end{eqnarray}

Using a Favre average is defined as
\begin{eqnarray}
   \mv{u}(\mv{x},t) = \mfa{\mv{u}}(\mv{x}) + \mv{\upp}
   \quad .
\end{eqnarray}


With the position vector \f$\mv{x} = (x,y,z) \f$. Applying the decomposition and the properties of the Reynolds operator (e.g. the mean of a perturbed value being zero \f$ \overline{u'}=0 \; \text{  or. } \; \mfa{u_i^{''}} = \overline{\rho u_i^{''}}=0  \f$ ) on the compressible Navier-Stokes equations one yields:

\f{align}{
\fracpart{\mta{\rho}}{t}+\fracpart{}{x_j}(\mta{\rho} \mfa{u_j}) &= 0 \\ 
\fracpart{\mta{\rho} \mfa{u_i}}{t}+\fracpart{}{x_j}(\mfa{u_j} \mta{\rho} \mfa{u_i}) &= -\fracpart{p}{x_i} + \fracpart{\overline{\sigma_{ij}}}{x_j} + \fracpart{\tau_{ij}}{x_j} \\ 
\fracpart{\mta{\rho} \mfa{E} }{t} + \fracpart{}{x_j}(\mfa{u_j}\mta{\rho} \mfa{H}) &=  \fracpart{}{x_j} (\overline{\sigma_{ij}} \mfa{u_i} + \overline{\sigma_{ij} u_i^{''}} ) - \fracpart{}{x_j} ( \mta{q_j} + c_p \overline{\rho u_j^{''} T''} - \mfa{u_i} \tau_{ij} + \frac{1}{2} \overline{\rho u_i^{''}  u_i^{''} u_j^{''}}) 
\f}

where  
\f{align}{
\mfa{H} &= \mfa{E} + \mta{p}/\mta{\rho} \\
\mta{q_j} &= -\overline{k_T \fracpart{T}{x_j}} \approx - \frac{c_p \mfa{\mu}}{Pr} \fracpart{\mfa{T}}{x_j} 
\f}

and the viscous stress tensor is defined as:
\f{align}{
\overline{\sigma_{ij}} &\approx 2\mfa{\mu} \left( \mfa{S_{ij}} - \frac{1}{3} \fracpart{\mfa{u_k}}{x_k} \delta_{ij} \right)
\f}

The derivation of these equations can (for example) be found in [[Gatski2009]]. 
The variable \f$ c_p\f$ describes the heat capacity at constant pressure, \f$ Pr = \frac{c_p \mu}{\kappa} \f$ is the Prandtl number with \f$ \kappa \f$ being the thermal conductivity. Sutherland's law is often used for computing the dynamic viscosity \f$ \mu\f$  
\f{align}{
 \mu = \mu_0 \left( \frac{T}{T_0} \right)^{3/2} \left( \frac{T_0 + C}{T+C} \right) ,
\f}

considering the relationship between the reference viscosity \f$ \mu_0 = 1.716 \times 10^{-5} \frac{kg}{m s} \f$, reference temperature \f$ T_0 = 273.15 K\f$, and the Sutherland constant \f$C = 110.4 K \f$. [[Sutherland1893]]  
The equation of state yields: 
\begin{eqnarray}
  \mta{p} = (\gamma-1) [\mta{\rho}\mfa{E} - \frac{1}{2} \mta{\rho} (\mfa{u}^2 +\mfa{v}^2 + \mfa{w}^2) - \mta{\rho}k ]
\end{eqnarray}
with the heat capacity ratio \f$\gamma \f$ and the turbulent kinetic energy \f$ k =\frac{1}{2}[(u_i^{''})^2 + (v_i^{''})^2+ (w_i^{''})^2 ]  \f$.  
Four terms emerge from the compressible RANS equations, that need to be modeled, namely the Reynolds stress term \f$ \tau_{ij} \f$, the turbulent heat flux \f$ c_p \overline{\rho u_j^{''} T''} \f$, the molecular diffusion term \f$ \overline{\sigma_{ij} u_i^{''}} \f$, and the turbulent transport term \f$ \frac{1}{2} \overline{\rho u_i'' u_i'' u_j'' } \f$.


