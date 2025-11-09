# Lagrangian particle tracking (LPT) # {#mmLPT}

## Motion Equation (Spherical Particles)
The equation of motion for a discrete particle \f$p\f$ with velocity \f$\mv{u}_p\f$ 
in non-dimensional form reads
\f{align}{ 
  \frac{d\mv{u}_p}{dt} &=
  (1 - \frac{\rho}{\rho_p} ) \frac{1}{Fr_0^2} \mv{e_g}
  + \frac{C_{D}}{24} \frac{Re_p}{\tau_p} (\mv{u} - \mv{u}_p), \\
  \frac{d\mv{x}_p}{dt} &= \mv{u}_p
\f}
with the Froude-number \f$Fr_0 = \frac{\mdim{a_{0}}}{\sqrt{\mdim{g}\, \mdim{L}}}\f$, 
unity gravity vector \f$\mv{e_g}\f$, particle Re number \f$Re_p\f$,
and particle relaxation time \f$\tau_p\f$ given as
\f{align}{
  Re_p = \frac{\rho||\mv{u} - \mv{u}_p|| d_p}{\mu}, \quad
  \tau_p = \frac{\rho_p d^2_p}{18\mu}.
\f}
The drag coefficient \f$C_D\f$ is given by a empirical relations, such as
\f{align}{
  C_D =
  \begin{cases}
    \frac{24}{Re'} & \text{for } Re'\leq 0.1\\
    \frac{24}{Re'}(1+\frac{1}{6} (Re') ^{\frac{2}{3}}) & \text{for } Re' \leq 1000\\
    0.424 & \text{for } Re' > 1000\\ 
  \end{cases}
\f}
from[[putnam1961integratable]] with Re number \f$Re' = Re_p Re_0\f$. 
The additional Reynolds number \f$Re_0\f$ appears
due to the non-dimensionalization of the fluid properties by 
the stagnation point reference value, i.e., \f$\mu = \frac{\mdim{\mu}}{\mdim{\mu_0}}\f$.

## Motion Equation (Ellipsoidal Particles)
The motion of a discrete rotational symmetric ellipsoidal particle \f$p\f$ with principal semi-minor axes
\f$a=b=c/\beta\f$ where \f$\beta\f$ refers to the aspect ratio is characterized by its velocity
\f$\mv{u}_p\f$ and angular velocity \f$\hat{\mv{\omega}}\f$.
Quantities in the particle fixed coordinate system are denoted by \f$\hat{\circ}\f$.

In the Stokes limit, the equation of motion for the translational movement of an ellipsoidal particle
[[Oberbeck1876]] reads
\f{align}{
  \frac{d\mv{u}_p}{dt} &=
  (1 - \frac{\rho}{\rho_p} ) \frac{1}{Fr_0^2} \mv{e_g}
  + \frac{24}{9} \frac{\beta^{2/3}}{\tau_p} \mm{R^{-1}}
  \begin{pmatrix}
  \frac{1}{\chi_0/a^2 + \alpha_0} & 0                               & 0 \\
  0                               & \frac{1}{\chi_0/a^2 + \alpha_0} & 0 \\
  0                               & 0                               & \frac{1}{\chi_0/a^2 + \beta^2 \gamma_0}
  \end{pmatrix}
  \mm{R} (\mv{u} - \mv{u}_p), \\
  \frac{d\mv{x}_p}{dt} &= \mv{u}_p
\f}
with the Froude-number \f$Fr_0 = \frac{\mdim{a_{0}}}{\sqrt{\mdim{g}\, \mdim{L}}}\f$,
unity gravity vector \f$\mv{e_g}\f$, rotation matrix \f$\mm{R}\f$, and the particle relaxation time \f$\tau_p\f$ given as
\f{align}{
  \tau_p = \frac{\rho_p d^2_p}{18 \rho \mu}.
\f}
\f$\alpha_0\f$, \f$\gamma_0\f$, and \f$\chi_0\f$ are shape factors that depend on the aspect ratio \f$\beta\f$.
Analytical expressions can be found in [[Happel&Brenner1965]].

For particles with \f$\beta=1\f$ the equation reduces to the simplified Maxey-Riley equation for spherical particles.

The rotational motion equation depends on the hydrodynamic torque [[Jeffrey1992]] and reads
\f{align}{
  \begin{pmatrix}
    \frac{d\omega_{\hat{x}}}{dt} \\
    \frac{d\omega_{\hat{y}}}{dt} \\
    \frac{d\omega_{\hat{z}}}{dt}
  \end{pmatrix} =
  \begin{pmatrix}
    \omega_{\hat{y}}\omega_{\hat{z}}\frac{\beta^2-1}{1+\beta^2} \\
    \omega_{\hat{z}}\omega_{\hat{x}}\frac{\beta^2-1}{1+\beta^2} \\
    0
  \end{pmatrix} +
  \frac{40}{9} \frac{\beta^{2/3}}{\tau_p}
  \begin{pmatrix}
    \frac{1}{\alpha_0 + \beta^2\gamma_0} \\
    \frac{1}{\alpha_0 + \beta^2\gamma_0} \\
    \frac{1}{2\alpha_0}
  \end{pmatrix}
  \begin{pmatrix}
    \frac{1-\beta^2}{1+\beta^2} \tau_{\hat{z}\hat{y}} + (\zeta_{\hat{z}\hat{y}}-\omega_{\hat{x}}) \\
    \frac{\beta^2-1}{1+\beta^2} \tau_{\hat{x}\hat{z}} + (\zeta_{\hat{x}\hat{z}}-\omega_{\hat{y}}) \\
    (\zeta_{\hat{y}\hat{x}}-\omega_{\hat{z}}) \\
  \end{pmatrix}
\f}
where the fluid shear stress and vorticity are defined as
\f{align}{
  \tau_{ij} = \frac{1}{2} \left( \fracpart{u_i}{x_j} + \fracpart{u_j}{x_i} \right), \quad
  \zeta_{ij} = \frac{1}{2} \left( \fracpart{u_i}{x_j} - \fracpart{u_j}{x_i} \right) .
\f}


## Energy Equation
For the total particle mass and temperature temporal rate of change 
(\f$\dot{m}_p\f$ ,\f$\dot{T}_p\f$ ) two
differential equations following the evaporation model of Miller 
and Bellan~([[miller1999]]-[[Miller1998]])
are included %in non-dimensional form
\f{align}{ \label{parteqm}
  \frac{d m_p}{d t} &= -\frac{Sh'}{3\,Sc' \, Re_0} (\frac{m_p}{\tau_p})\ln(1+B_{M,neq})\\
  \label{parteqt}
    \frac{dT_p}{dt} &=\frac{Nu'_{neq}}{3\,Pr'\,Re_0}\theta_1 \frac{T-T_p}{\tau_p}
    +L_{ev} (\gamma - 1) \theta_1 \frac{\dot{m_p}}{m_p}
        \quad ,
\f} 
where \f$\theta_1\f$ is the ratio of fluid-to-liquid heat capacity 
\f$\theta_1=\frac{\mdim{C_p}}{\mdim{C_{p,p}}}\f$ and
\f$L_{ev}=\frac{\mdim{L_{ev}}}{\mdim{a_{0}}}\f$ the non-dimensional latent heat of vaporization.
The Schmidt and Prandtl number of the fluid phase are given by
\f{align}{ \label{eqSchmidt}
Sc' &= \frac{\mu}{\rho \Gamma} Sc_0 &\quad \text{with} \quad 
Sc_0 =\frac{\mdim{\mu_0}}{\mdim{\rho_0} \mdim{\Gamma_0}} \\
Pr' &= \frac{\mu C_p}{\lambda} Pr_0 &\quad \text{with} \quad 
Pr_0 = \frac{\mdim{\mu_0} \mdim{C_p}}{\mdim{\lambda_0}} \\
\f} 
with the dimensionless binary diffusion coefficient \f$\Gamma\f$, the
thermal conductivity \f$\lambda\f$, and the respective dimensional reference state values.
For the convective heat and mass transfer, the empirical Ranz-Marshall correlations
([[RanzMarshall1]]-[[RanzMarshall2]])
for the Sherwood and Nusselt number
\f{align}{ \label{eqRanzMarshall}
 Sh' &= 2+0.552 Re'^{\frac{1}{2}} Sc'^{\frac{1}{3}} \\
 Nu' &= 2+0.552 Re'^{\frac{1}{2}} Pr'^{\frac{1}{3}}
\f} 
are used. 
The equilibrium vapor mole fraction is obtained from the Clausius-Clapeyron equation
\f{align}{
  \chi_{s,eq} =
  \frac{p_{B}}{p}\exp\left(\frac{\gamma L_{ev}}{\theta_2}\left(\frac{1}{T_B}-\frac{1}{T_p}\right) \right)
  \label{eq::clausius}
\f}
with the liquid boiling point (\f$T_B\f$ at \f$p_B\f$) and the mole-fraction ratio \f$\theta_2 = \frac{M}{M_p}\f$.
Non-equilibrium effects are considered by a correction of 
the Nusselt number with the non-dimensional evaporation parameter \f$\beta\f$ 
and by a correction of the vapor mole fraction at the surface by the Langmuir-Knudsen law
\f{align}{ \label{eqNeq}
  Nu'_{neq} &= Nu' \frac{\beta}{e^{\beta}-1} &\quad{\text{with} } \quad  &\beta = 0.5 Pr' Re_b'\\
  \chi_{s,neq} &= \chi_s - 2 \beta L_k, &\quad{ \text{with} } \quad 
  &L_k = \frac{\mu \sqrt{2 \pi T_p \frac{\theta_2}{\gamma}}}{p d_p Sc' Re_0}
\f} 
with \f$Re_b' = \frac{\dot{m}_p}{\mu \pi d_p} Re_0\f$ as  
the Reynolds number based on the blowing velocity \f$u_b\f$, 
which is obtained from mass conservation at the droplet surface 
with \f$\dot{m}_p=-\pi d_p^2 u_b \rho\f$ as in[[miller1999]].

The non-equilibrium species vapor mass fraction at the particle surface \f$Y_{s,neq}\f$ and the
Spalding transfer number for mass \f$B_{M,neq}\f$ are defined as
\f{align}{ \label{eqNeq2}
Y_{s,neq} =  \frac{\chi_{s,neq}}{\chi_{s,neq} + (1-\chi_{s,neq})\theta_2} \, \text{and} \quad  
B_{M,neq} = \frac{Y_{s,neq} - Y}{1 - Y_{s,neq}}. 
\f} 
The "1/3 rule" following [[Hubbard]] is used as a reference state 
in the boundary layer around the particle to obtain fluid properties. 
The Wilke mixing rule [[Wilke]] is applied for the fuel vapor and 
fluid mixing in the particle boundary layer.

The following source terms are used for the mass, momentum and 
energy exchange between the particle and the fluid phase
\f{align}{
  \begin{bmatrix}
    S_\rho \\
    \mv{S_u} \\[6pt]
    S_E \\[6pt]
    S_Y
  \end{bmatrix} = 
  \begin{bmatrix}
  \dot{m}_p \\
  m_p \frac{C_{D}}{24} \frac{Re_p}{\tau_p Re_0} (\mv{u} - \mv{u}_p) 
   + \dot{m}_p \mv{u}_p  \\[6pt]
  \mv{S_i} \cdot \mv{u}_p - 0.5 \dot{m}_p \mv{u}_p \cdot \mv{u}_p 
  + \frac{C_{p,p}}{\gamma -1} \left(m_p \dot{T}_p + \dot{m}_p T_p \right) \\[6pt]
  \dot{m}_p
  \end{bmatrix}.
\f}

## Spray Modeling
To reduce the computational effort generated by the spray injection, 
a number of fuel droplets \f$N_d\f$ are grouped into a parcel. 
All droplets in a parcel are assumed to have the same physical values and 
suffer from the same change in their properties, i.e., they are tracked 
as a single Lagrangian particle.

The spray atomization and secondary breakup process, 
leading to decreasing particle diameter, is modeled by a combination 
of the Kelvin-Helmholtz (KH) and Rayleigh-Taylor (RT) model.
The present implementation is based on the formulation of 
Reitz [[reitz1987]] and Beale and Reitz [[Reitz1999]].

### Kelvin Helmholtz Model
The Kelvin-Helmholtz (KH) wave model assumes droplet
break up due to wave growth excited by aerodynamic forces based 
on the relative velocity between liquid and gas phase. 
Continuous mass stripping from the parent droplet into 
smaller child droplets occurs.
A correlation for the most unstable wave growth rate 
\f$\Omega_{KH}\f$ and  wavelength \f$\Lambda_{KH}\f$ can be 
obtained from first order linear analysis such that
\f{align}{
  \Omega_{KH} &=
  \frac{0.34 + 0.38 We'^{1.5}}{(1 + Z')(1 + 1.4 T'^{0.6})}
  \sqrt{\frac{8 \sigma_p}{\rho_p d^3_p We_0}} \\
  \Lambda_{KH} &= \frac{4.51 d_p (1 + 0.45 \sqrt{Z'})(1 + 0.4 T'^{0.7})}
  {(1 + 0.865 We'^{1.67})^{0.6}},
\f} 
with the Ohnesorge number \f$Z'={\sqrt{We_l'}}/{Re_{l}'}\f$ and
the Taylor number \f$T' = Z'\sqrt{We'}\f$. 
The Weber numbers are defined as
\f$We'=\frac{0.5 \rho ||\mv{u}-\mv{u}_p||^2 d_p}{\sigma_p} We_0\f$, 
the liquid Weber number is
\f$We_l'=\frac{0.5 \rho_p||\mv{u}-\mv{u}_p||_2^2d_p}{\sigma_p} We_0\f$
and the liquid Reynolds number is
\f$Re_{l}' = 0.5 \frac{\rho_p||\mv{u} - \mv{u}_p||_2 d_p}{\mu_p} Re_0\f$. 
The breakup time \f$\tau_{KH}\f$ and 
droplet size \f$d_{KH}\f$ can be determined as
\f{align}{
  \tau_{KH} &= \frac{1.863 B_1 d_p}{\Omega_{KH} \Lambda_{KH}} \quad \text{and} \quad
  d_{KH} = 2 B_0 \Lambda_{KH}.
\f}
The model adjustable control parameters are \f$B_0\f$ and \f$B_1\f$.
If a KH breakup occurs, the shedded mass of child droplets is
accumulated until a mass limit \f$m_{KH,lim}\f$
or the parent droplet mass has been reached.
If the mass limit is reached, 
a new parcel with the mass of the shedded droplets and \f$d_{KH}\f$
is added. 
The number of droplets \f$N_d\f$ in the parent parcel 
remains the same, and a new parent diameter is computed 
based on mass conservation.
If the parent mass is reached, the parent diameter is set to 
\f$d_{KH}\f$ and the number of droplets is adapted such that mass conservation is achieved.
See Aguerre et. al.[[Aguerre]] for further details of this implementation.

### Rayleigh Taylor Model
The Rayleigh-Taylor (RT) model describes the catastrophic breakup mechanism
caused by deceleration effects on the liquid droplet due to aerodynamic drag.
The growth rate \f$\Omega_{RT}\f$ and wavelength \f$\Lambda_{RT}\f$ of the RT
wave are defined as
\f{align}{
  \Omega_{RT} &= \sqrt{\frac{2}{3 \sqrt{3\sigma_p}}
  \frac{(||\mv{a}_p|| (\rho_p - \rho))^{\frac{3}{2}}}{\rho_p + \rho} } We_0^{0.25}\\
  \Lambda_{RT} &= 2 \pi \sqrt{\frac{3 \sigma_p}
        {||\mv{a}_p|| (\rho_p - \rho)}},
\f} 
with \f$\mv{a}_p\f$ being the droplet acceleration and 
\f$\sigma_p\f$ the surface tension. 
The droplet size after breakup \f$d_{RT}\f$ and 
the RT breakup time \f$\tau_{RT}\f$ are given by
\f{align}{
  \tau_{RT} =\frac{1}{\Omega_{RT}} \quad \text{and} \quad
  d_{RT} = C_{RT} \Lambda_{RT} \label{eq_RT}
\f} 
with \f$C_{RT}\f$ as the primary model control parameter.
If RT breakup occurs the parcel diameter is reduced to \f$d_{RT}\f$ and 
the number of droplets is increased based on mass conservation.

The interaction of the two breakup models is based on 
Beale and Reitz[[Reitz1999]],
where it is assumed that within a certain distance \f$L_{bu}\f$ 
the spray consists of a liquid core
\f{align}{
L_{bu} = \frac{C_{bu}}{2} d_i \sqrt{\frac{\rho_p}{\rho_i}}
\quad ,
\f}
where \f$C_{bu}\f$ is a model parameter, 
\f$d_i\f$ the injector orifice diameter, and 
\f$\rho_i\f$ the initial fluid density.
Within \f$L_{bu}\f$ only KH breakup is considered, while beyond 
both mechanisms are applied 
with the RT model as the predominant mechanism.

## References


# Empirical Drag law correlations # {#mmLPTDragLaw}


