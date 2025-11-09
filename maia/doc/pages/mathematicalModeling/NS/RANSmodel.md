# RANS model  # {#mmRANSmodel}

## Modeling approach 
Their exist a wide range of models approximating the unknown terms. Most commonly a Reynolds analogy is chosen
for modeling the turbulent heat flux, with the turbulent Prandtl number \f$Pr_t=0.9\f$ and the assumptions of
an ideal gas for determining the heat capacity at constant pressure \f$c_p \f$.

\f{align}{
  c_p \overline{\rho u_j'' T''} \approx - \frac{c_p \mfa{\mu_t}}{Pr_t} \fracpart{\mfa{T}}{x_j}
  \f}

The molecular diffusion and turbulent transport expressions can be lumped together and might be modeled by:
\f{align}{
\overline{\sigma_{ij} u_i^{''}} - \frac{1}{2} \overline{\rho u_i'' u_i'' u_j'' } \approx \left( \mfa{\mu} + \frac{\mfa{\mu_t}}{\sigma_k} \right) \fracpart{k}{x_j}
  \f}
where the coefficient \f$ \sigma_k\f$ is associated with the modeling equation for k (and is sometimes neglected).  

<!--The mathematical background of the ansatz chosen for m-AIA will be presented briefly in the following. --> 
The Boussinesq approximation is commonly used to model the Reynolds stress:
\f{align}{
  \tau_{ij} = 2\mfa{\mu_t} \left( \mfa{S_{ij}} -\frac{1}{3} \fracpart{\mfa{u_k}}{x_k} \delta_{ij} \right) - \frac{2}{3} \mta{\rho}k \delta_{ij}
\f}
The last term is ignored as k is not readily available. The eddy viscosity \f$\mfa{\mu_t} \f$ as well as the term
\f$\mfa{S_{ij}} = \frac{1}{2}(\fracpart{\mfa{u_i}}{x_j} + \fracpart{\mfa{u_j}}{x_i}) \f$ are determined
from the turbulence model, for m-AIA a compressible form of the Spalart-Allmaras One-Equation Model (SA-noft2-Catris) [[Spalart-Allmaras]].



###A Compressible Form of Spalart-Allmaras One-Equation Model

The one-equation model in compressible form is given by the following equation: 

\f{align}{ \frac{\partial \rho \tilde \nu}{\partial t} + u_j \frac{\partial \rho \tilde \nu}{\partial x_j} = c_{b1} \tilde S \rho \tilde \nu - c_{w1} f_{\omega} \rho \left(\frac{\tilde \nu}{d} \right)^2 + \frac{1}{\sigma} \frac{\partial}{\partial x_j} \left(\mu \frac{\partial \tilde \nu}{\partial x_j} \right) + \frac{1}{\sigma} \frac{\partial}{\partial x_j} \left(\sqrt \rho \tilde \nu \frac{\partial \sqrt \rho \tilde \nu}{\partial x_j} \right) + \frac{c_{b2}}{\sigma} \frac{\partial \sqrt \rho \tilde \nu}{\partial x_i} \frac{\partial \sqrt \rho \tilde \nu}{\partial x_i}\f}

The turbulent eddy viscosity is computed by

\f{align}{ \hat \mu_t = \rho \tilde \nu f_{v1} 
\quad . \f}

where

\f{align}{ f_{v1} = \frac{\chi^3}{\chi^3 + c_{v1}^3} \f}

\f{align}{ \chi = \frac{\tilde \nu}{\nu}  \f}

\f{align}{ \tilde S = \Omega + \frac{\tilde \nu}{\kappa^2 d^2} f_{v2}  \f}

\f{align}{ \Omega = \sqrt{2 W_{ij} W_{ij}} \f}

\f{align}{ f_{v2} = 1 - \frac{\chi}{1 + \chi f_{v1}}  \f}

\f{align}{ f_{\omega} = g \left[\frac{1 + c_{w3}^6}{g^6 + c_{w3}^6}\right]^{1/6}  \f}

\f{align}{ g = r + c_{w2} \left(r^6 -r \right)  \f}

\f{align}{ r = min \left[\frac{\tilde \nu}{\tilde S \kappa^2 d^2}, 10 \right] \f}

\f{align}{ f_{t2} = c_{t3} exp \left(-c_{t4} \chi^2 \right) \f}

\f{align}{ W_{ij} = \frac{1}{2} \left(\frac{\partial \hat u_i}{\partial x_j} - \frac{\partial \hat u_j}{\partial x_i} \right) \quad . \f} 

The constants are:

\f$ c_{b1} = 0.1355 \f$ 
\f$ \sigma = 2/3 \f$    
\f$ c_{b2} = 0.622 \f$  
\f$ \kappa = 0.41 \f$   
\f$ c_{w2} = 0.3 \f$    
\f$ c_{w3} = 2 \f$      
\f$ c_{v1} = 7.1 \f$    
\f$ c_{t3} = 1.2 \f$    
\f$ c_{t4} = 0.5 \f$    
\f$ c_{w1} = \frac{c_{b1}}{\kappa^2} + \frac{1 + c_{b2}}{\sigma} \f$    
