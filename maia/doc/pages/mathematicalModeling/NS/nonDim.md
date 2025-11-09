# Nondimensionalization of the Navier-Stokes Equations # {#mmNonDimensionalization}

Within @maia the Navier-Stokes equations are solved using
nondimensionalized variables such that the numerical solution is
valid for all fluid flows with dynamic similarity, i.e., with identical
similarity parameters. The nondimensional form is obtained by
normalizing all dimensional variables denoted by the superscript \f$ ^* \f$ with reference values denoted by the index \f$ _{ref} \f$. The reference
values used in the @maia formulation of the Navier-Stokes equations are based on stagnation point values, for which an isentropic deceleration of the fluid to a velocity of zero is assumed and which are denoted by an index \f$ _0 \f$. In the table below the definition of nondimensional values is summarized.

Quantity             | nondimensional variable |  Reference Value 
---------------------|-------------------------|------------------------  
length               |\f$ \mv{x}=\frac{\mv{x}^*}{L^*_{ref}} \f$ |\f$ L^*_{ref} = L^* \f$         
velocity             |\f$ \mv{u}=\frac{\mv{u}^*}{u^*_{ref}} \f$ |\f$ u^*_{ref} = a^*_0\f$        
time                 |\f$ t=\frac{t^*}{t^*_{ref}} \f$	        |\f$ t^*_{ref} = \frac{L^*}{a^*_0}\f$        
density              |\f$ \rho=\frac{\rho^*}{\rho^*_{ref}} \f$|\f$ \rho^*_{ref} = \rho^*_0 \f$  
pressure             |\f$ p=\frac{p^*}{p^*_{ref}} \f$	        |\f$ p^*_{ref} = \rho^*_0 {a^*_0}^2\f$        
temperature          |\f$ T=\frac{T^*}{T^*_{ref}} \f$		|\f$ T^*_{ref} = T^*_0\f$        
energy               |\f$ e=\frac{e^*}{{u^*_{ref}}^2} \f$	|\f$ {u^*_{ref}}^2 = {a^*_0}^2 \f$
speed of sound       |\f$ a_{ref} \f$                           |\f$ a_0 \f$
thermal conductivity |\f$ k=\frac{k^*}k^*_{ref} \f$  	        |\f$ k^*_{ref} = k^*_0\f$        
viscosity            |\f$ \mu=\frac{{\mu}^*}{\mu^*_{ref}} \f$   |\f$ \mu^*_{ref} =\mu^*_0 \f$      


The reference length \f$ L^* \f$ is defined by the geometry data provided in the STL geometry file and corresponds to the length, where the difference of any two points of the coordinates \f$ \Delta \mv{x} \f$ becomes unity.

With the above choice of reference values, the following nondimensional similarity parameters are obtained:

Similarity parameter | Value
---------------------|------
Strouhal number      | \f$ Sr = \frac{L^*_{ref}}{u^*_{ref} t^*_{ref}} = 1 \f$  
Reynolds number      | \f$ Re = \frac{\rho^*_{ref} u^*_{ref} L^*_{ref}}{\mu^*_{ref}} = \frac{\rho^*_0 a^*_0 L^*}{\mu^*_0} \f$  
Euler number         | \f$ Eu = \frac{p^*_{ref}}{\rho^*_{ref} {u^*_{ref}}^2} = 1  \f$ 
Prandtl number       | \f$ Pr = \frac{\mu^*_{ref} {c^*_p}_{ref}}{k^*_{ref}}\f$   
Mach number          | \f$ M = \frac{u^*_{ref}}{a^*_{ref}} = 1  \f$ 


For this particular choice of reference, i.e., stagnation point
values, the Strouhal, the Euler, and the Mach number attain a value of
unity such that they do not appear in the nondimensional form of the
Navier-Stokes equations in @maia. The other similarity parameters, i.e., the
Reynolds number, the Prandtl number, and the ratio of specific heats are used in various functions for the computation of convective and viscous fluxes **add links**.

For the determination of the pressure and temperature, the fluid is assumed
to be an **ideal gas** with constant specific heats. Then, the speed
of sound \f$ a \f$ can be computed as

\begin{equation}
  a^*_{ref} = \sqrt{\gamma R^* T^*_{ref}} = a^*_0 = \sqrt{\gamma R^* T^*_0}
\end{equation}

with the ideal gas constant \f$ R^* \f$ and the isentropic coefficient
or ratio of specific heats \f$ \gamma = \frac{c^*_p}{c^*_v} \f$. The
specific heat can be computed as, e.g.,

\begin{equation}
  c^*_p = \frac{\gamma}{\gamma -1} R^*
  \quad .
\end{equation}

The nondimensional speed of sound \f$ a \f$ is then defined by

\begin{equation}
  a = \sqrt{T}
  \quad .
\end{equation}

The nondimensional temperature \f$ T \f$ can be determined by the ideal gas law in nondimensional form

\begin{equation}
  T = \gamma \frac{p}{\rho}
  \quad .
\end{equation}

For an ideal gas, the nondimensional pressure \f$ p \f$ can be computed as:
\begin{equation}
  p = (\gamma - 1 )(\rho E - \frac12 \rho \mv{u}^2 )
\end{equation}


Note that the nondimensional pressure at a stagnation point \f$ p_0 \f$ for an isentropic deceleration to a velocity of zero becomes

\begin{equation}
  p_0 = \frac{1}{\gamma}
  \quad .
\end{equation}

Likewise, the nondimensional total temperature \f$ T_0 \f$ becomes

\begin{equation}
  T_0 = 1
  \quad .
\end{equation}


For external flows, often reference values based on undisturbed free
stream values, denoted by the index \f$ _\infty \f$, are used. In the
following, the relations between stagnation point and free stream
values are summarized. If the Mach number \f$ M_{\infty} \f$ and the
Reynolds number \f$ Re_{\infty} \f$ of the free stream flow are given,
then the Reynolds number based on stagnation point values can be
computed as follows

\begin{equation}
  Re_0 =  \frac{\rho^*_0 a^*_0}{\mu^*_0} \frac{\mu^*_{\infty}}{\rho^*_{\infty} u^*_{\infty}} Re_{\infty}
  \quad . 
\end{equation}

With the relations for **isentropic flows** of an **ideal gas**

\begin{equation}
\frac{T^*_{\infty}}{T^*_0} = \left( 1 + \frac{\gamma - 1}{2} M_{\infty}^2 \right)^{-1}
\end{equation}

\begin{equation}
\frac{p^*_{\infty}}{p^*_0} = \left( \frac{\rho^*_{\infty}}{\rho^*_0} \right)^{\gamma} = \left( \frac{T^*_{\infty}}{T^*_0} \right)^{\frac{\gamma}{\gamma -1}}
\end{equation}

the following two expressions are obtained
\begin{equation}
\frac{\rho^*_0}{\rho^*_{\infty}} = \left( \frac{T^*_0}{T^*_{\infty}} \right)^{\frac{1}{\gamma - 1}} = \left( 1 + \frac{\gamma - 1}{2} M_{\infty}^2 \right)^{\frac{1}{\gamma - 1}}
\end{equation}


\begin{equation}
\frac{a^*_0}{u^*_{\infty}} = \frac{a^*_{\infty}}{u^*_{\infty}} \frac{a^*_0}{a^*_{\infty}} =
\frac{1}{M_{\infty}} \sqrt{\frac{T^*_0}{T^*_{\infty}}} =
\frac{1}{M_{\infty}} \left( 1 + \frac{\gamma - 1}{2} M_{\infty}^2 \right)^{\frac{1}{2}}
\quad .
\end{equation}

When a power law for the viscosity is assumed for simplicity

\begin{equation}
\frac{\mu^*_{\infty}}{\mu^*_0} = \left( \frac{T^*_{\infty}}{T^*_0} \right)^{\omega} = \left( 1 + \frac{\gamma - 1}{2} M_{\infty}^2 \right)^{-\omega}
\quad ,
\end{equation}

where \f$ \omega \f$ for air is approximately \f$ \omega \approx 0.72 \f$, one finally obtains

\begin{equation}
  Re_0 = \frac{Re_{\infty}}{M_{\infty}}  \left( 1 + \frac{\gamma - 1}{2} M_{\infty}^2 \right)^{\frac12 \frac{\gamma+1}{\gamma-1} - \omega}
  \quad .
\end{equation}

For air, i.e., \f$ \gamma = 1.4 \f$, and \f$ \omega = 0.72 \f$, the following factors for the Reynolds number conversion are obtained:

free stream Mach number \f$ M_{\infty} \f$ | Reynolds number \f$ Re_0 \f$
-------------------------------------------|----------------------------
\f$ M_{\infty} =0.01 \f$        | \f$ Re_0 = 100.06 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.1 \f$        | \f$ Re_0 = 10.58 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.2 \f$        | \f$ Re_0 = 6.21 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.3 \f$        | \f$ Re_0 = 5.30 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.4 \f$        | \f$ Re_0 = 5.38 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.5 \f$        | \f$ Re_0 = 6.05 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.6 \f$        | \f$ Re_0 = 7.20 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.7 \f$        | \f$ Re_0 = 8.85 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.8 \f$        | \f$ Re_0 = 11.04 \, Re_{\infty} \f$ 
\f$ M_{\infty} =0.9 \f$        | \f$ Re_0 = 13.86 \, Re_{\infty} \f$ 
\f$ M_{\infty} =1.0 \f$        | \f$ Re_0 = 17.40 \, Re_{\infty} \f$ 
\f$ M_{\infty} =2.0 \f$        | \f$ Re_0 = 118.40 \, Re_{\infty} \f$ 
\f$ M_{\infty} =3.0 \f$        | \f$ Re_0 = 445.57 \, Re_{\infty} \f$ 


## Additional quantities

The dimensional equation for the wall-shear stress 
\f{align}{
 \tau^*_w = \mu^* \frac{\partial u^*}{\partial y^*}
\f}

is nondimensionalized as

\f{align}{
\tau_w = \frac{\tau_w^*}{a_0^{*2} \rho^*_0}
\f}

such that

\f{align}{
  \tau_w = \frac{1}{a_0^{*2} \rho^*_0} \frac{\mu^*_0 a^*_0}{L^*} \mu \frac{\partial u}{\partial y} =  \underbrace{\frac{\mu^*_0}{a^*_0 \rho^*_0 L^*}}_{= 1/Re_0} \mu \frac{\partial u}{\partial y} = \frac{1}{Re_0} \mu \frac{\partial u}{\partial y}
  \quad .
\f} 

For instance, the skin-friction coefficient \f$ c_f \f$ is computed from nondimensionalized values by

\f{align}{
  c_f  = \frac{\tau_w}{\frac{\rho_\infty}{2}u_\infty^2} = \frac{\mu \frac{\partial u}{\partial y}}{Re_0 \frac{\rho_\infty}{2}u_\infty^2}
  \quad .
\f}

The dimensional equation for the friction velocity \f$ u^*_\tau = \sqrt{\tau^*_w / \rho^*}\f$ is nondimensionalized as

\f{align}{
   u_\tau = \frac{u^*_\tau}{a^*_0} = \frac{1}{a^*_0} \sqrt{\frac{\mu \frac{\partial u}{\partial y} a_0^{*2} \rho^*_0}{Re_0 \rho \rho^*_0}} = \sqrt{\frac{\mu \frac{\partial u}{\partial y} }{Re_0 \rho }} = \sqrt{\frac{\mu \frac{\partial u}{\partial y} }{Re_0 \rho}}
   \quad .
\f}

Inner scaling (also known as viscous scaling), i.e., scaled with the friction velocity \f$u_\tau\f$ and the kinematic viscosity \f$\nu\f$ can be derived also for nondimensional quantities. In the following the nondimensionalization is shown for several quantities.

* Velocities
        \f{align}{
        u^+ = \frac{u^*}{u^*_\tau} = \frac{u a^*_0}{u_\tau a^*_0} = \frac{u}{u_\tau}
        \f}
* Lengths
 \f{align}{
  y^+  = \frac{y^* u^*_\tau \rho^*}{\mu^*} = \frac{y L^* u_\tau a^*_0 \rho \rho^*_0}{\mu \mu^*_0} = \frac{y~ u_\tau \rho}{\mu}  \underbrace{\frac{L^* a^*_0 \rho^*_0}{\mu^*_0}}_{Re_0} = Re_0 \frac{y~ u_\tau \rho}{\mu}
  \f}
* Time
\f{align}{
  t^+  = \frac{t^* u^{*2}_\tau \rho^*}{\mu^*} = \frac{t \frac{L^*}{a^*_0} u^2_\tau a^{*2}_0 \rho \rho^*_0}{\mu \mu^*_0} = \frac{t~ u^2_\tau \rho}{\mu}  \underbrace{\frac{L^* a^*_0 \rho^*_0}{\mu^*_0}}_{Re_0} = Re_0 \frac{t~ u^2_\tau \rho}{\mu}
  \f}
* Vorticity
 \f{align}{
 \omega^+ = \frac{\mu^* \omega^*}{\rho^* u^{*2}_\tau} = \frac{\mu \omega \mu^*_0 a^*_0}{\rho u^2_\tau \rho^*_0 a^{*2}_0 L^*} = \frac{\mu \omega}{u^2_\tau Re_0}
\f}



