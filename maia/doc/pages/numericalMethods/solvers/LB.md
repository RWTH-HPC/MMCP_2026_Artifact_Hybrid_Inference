# Lattice Boltzmann (LB) # {#nmLB}

[TOC]

## Lattice Boltzmann equation
Discretizing the Boltzmann equation according to [Chen1992] and [Qian1992] yields the lattice Boltzmann equation
\f{equation}{
  f_i(\mv{x} + \mv{\xi}_i \Delta t, t + \Delta t) = f_i(\mv{x}, t) + \omega (f_i^{eq}(\mv{x}, t)-f_i(\mv{x}, t)) \,,
  \label{eq:lattice_boltzmann}
\f}
with the discrete PPDF \f$ f_i \f$, the time increment \f$ \Delta t \f$ and the discrete particle velocity \f$ \mv{\xi}_i \f$ in the discrete direction \f$ i\f$.
For time integration, a simple Euler forward scheme is used.

The discrete Maxwell distribution function is defined by
\f{equation}{
  f_i^{eq} = t_p \rho \left[ 1 + \frac{\mv{\xi}\cdot\mv{u}}{c_s^2}
  + \frac{1}{2} \left( \frac{\mv{\xi} \cdot \mv{u}}{2c_s^2} \right)^2 - \frac{\mv{u}^2}{c_s^2} \right]
  \label{eq:equilibrium}
  \,,
\f}

where \f$ c_s=1/\sqrt{3} \f$ is the isothermal speed of sound.
The coefficient \f$ t_p \f$ is given for the most popular discretization models mentioned by [Qian1992] in the table below. 

<table>
<tr><th> @f$ t_p @f$ <th> D2Q9 <th> D3Q19 <th> D3Q27
<tr><td> @f$ \frac{\mv{\xi}_i^2}{\mv{\xi}_0^2}=0 @f$ <td> @f$ \frac{4}{9} @f$ <td> @f$ \frac{1}{3} @f$ <td> @f$ \frac{8}{27} @f$
<tr><td> @f$ \frac{\mv{\xi}_i^2}{\mv{\xi}_0^2}=1 @f$ <td> @f$ \frac{1}{9} @f$ <td> @f$ \frac{1}{18} @f$ <td> @f$ \frac{2}{27} @f$
<tr><td> @f$ \frac{\mv{\xi}_i^2}{\mv{\xi}_0^2}=2 @f$ <td> @f$ \frac{1}{36} @f$ <td> @f$ \frac{1}{36} @f$ <td> @f$ \frac{1}{54} @f$
<tr><td> @f$ \frac{\mv{\xi}_i^2}{\mv{\xi}_0^2}=3 @f$ <td> - <td> - <td> @f$ \frac{1}{216} @f$
</table>

The denotation DxQy describes the spatial dimension (x) and the number of molecular velocities (y).
Molecules are only allowed to move from a node to a neighbouring node (\f$ \Delta x \f$) 
in on time step (\f$ \Delta t \f$), the mean molecular velocity can also be defined as \f$ \xi_0=\Delta x / \Delta t \f$. 

The algorithm to solve Eq.\f$\eqref{eq:lattice_boltzmann}\f$ cosists of two steps:

-   Collision step:

@f$\tilde{f_i}\left(\mv{x},t\right) = f_i\left(\mv{x},t\right)+ \omega \left(f_i^{eq}\left(\mv{x},t\right)-f_i\left(\mv{x},t\right)\right)@f$

-   Propagation step:

@f$f_i\left(\mv{x}+\xi_i\Delta t,t+\Delta t\right)=\tilde{f_i}\left(\mv{x},t\right)@f$

Using the Chapman-Enskog expansion, the relation between \f$ \omega \f$ and the kinematic shear viscosity \f$ \nu \f$ can be derived as
\f{equation}{
  \nu = c_s^2 \left(\frac{\Delta t}{\omega}-\frac{\Delta t}{2}\right).
  \label{eq:omega_nu}
\f}
Considering the limits of Eq. \f$\eqref{eq:omega_nu}\f$, 
the stability limits of the LB method are shown. 
Since the viscosity only has a physical meaning for \f$ \nu>0 \f$, 
the upper limit of the relaxation frequency is \f$ \omega<2  \f$. 
For a small collision frequency, i.e., \f$ \omega \rightarrow 0 \f$, 
Eq. \f$\eqref{eq:omega_nu}\f$ yields a large viscosity \f$ \nu \rightarrow \inf \f$.
If \f$ \omega \rightarrow 2 \f$, the LB method tends to become unstable.

## Thermal lattice Boltzmann method
The LB scheme described above is valid for low Mach numbers,
and hence the energy equation decouples from the momentum
equations. To simulate the temperature distribution, additional
sets of PPDFs are solved on the same lattice.
Three approaches are implemented for computing the temperature.
The basic thermal approach is based on [Shan1997] and the temperature distribution function \f$ g_i \f$ is calculated by
\f{equation}{
  g_i(\mv{x} + \mv{\xi}_i \Delta t, t + \Delta t) = g_i(\mv{x}, t) + \omega_g (g_i^{eq}(\mv{x}, t)-g_i(\mv{x}, t)) \,.
  \label{eq:lattice_boltzmann_thermal}
\f}

The thermal relaxation frequency \f$ \omega_g \f$ is a function of the thermal diffusivity \f$ \kappa \f$
\f{equation}{
  \kappa = c_s^2 \left(\frac{\Delta t}{\omega_g}-\frac{\Delta t}{2} \right) \,.
  \label{eq:omega_g}
\f}

The thermal equilibrium distribution function \f$ g_i^{eq} \f$ for the basic approach is given by
\f{equation}{
  g_i^{eq} = t_p T \left[ 1 + \frac{\mv{\xi}\cdot\mv{u}}{c_s^2}
  + \frac{1}{2} \left( \frac{\mv{\xi} \cdot \mv{u}}{2c_s^2} \right)^2 - \frac{\mv{u}^2}{c_s^2} \right]
  \label{eq:equilibrium_thermal}
  \,,
\f}

With their internal energy LB method [He1998] improved the method of [Shan1997]. 
That is, the internal energy \f$ \epsilon = DRT/2 \f$ multiplied by \f$ \rho \f$ is transported by the internal energy distribution function \f$ h_i \f$
\f{equation}{
  h_i(\mv{x} + \mv{\xi}_i \Delta t, t + \Delta t) = h_i(\mv{x}, t) + \omega_h (h_i^{eq}(\mv{x}, t)-h_i(\mv{x}, t)) \,.
  \label{eq:lattice_boltzmann_inner}
\f}

The relation between the thermal relaxation frequency \f$ \omega_h \f$ to \f$ \kappa \f$ differs from Eq. \f$\eqref{eq:omega_g}\f$
\f{equation}{
  \kappa = c_s^2 \frac{D+2}{D} \left(\frac{\Delta t}{\omega_h}-\frac{\Delta t}{2} \right) \,.
  \label{eq:omega_h}
\f}

The thermal equilibrium reads
\f{equation}{
  h_i^{eq} = t_p \rho \epsilon \left[ \frac{\mv{\xi}^2}{Dc_s^2} 
  + \left(\frac{\mv{\xi}^2}{Dc_s^2} - \frac{2}{D} \right)
  \times \frac{\mv{u} \cdot \mv{\xi}}{c_s^2} + \frac{\mv{u}\cdot\mv{u}}{2c_s^2}\cdot\left(\frac{\mv{\xi} \cdot \mv{\xi}}{c_s^2}-\delta\right)\right] 
  \label{eq:equilibrium_inner}
  \,,
\f}
where \f$ \delta \f$ stands for the Kronecker delta.

The third thermal approach considers the total energy \f$ E = e + 1/2 u^2 \f$ according to [Guo2007], 
where \f$ u \f$ is the velocity magnitude.
The total energy distribution function is formulated as
\f{equation}{
  k_i(\mv{x} + \mv{\xi}_i \Delta t, t + \Delta t) = k_i(\mv{x}, t) + \omega_k (k_i^{eq}(\mv{x}, t)-k_i(\mv{x}, t)) + \\
  (\omega_k - \omega) \left(\mv{\xi} \cdot \mv{u} - \frac{\mv{u} \cdot \mv{u}}{2} \right) (f_i(\mv{x}, t)-f_i^{eq}(\mv{x}, t)) \,.
  \label{eq:lattice_boltzmann_total}
\f}

The collision frequency is in this case related to the thermal conductivity by
\f{equation}{
  \kappa = c_s^2 \left(\frac{\Delta t}{\omega_k}-\frac{\Delta t}{2} \right) \,.
  \label{eq:omega_k}
\f}

The equilibrium distributions are given by
\f{equation}{
  k_i^{eq} = t_p \rho c_s^2 \left[ \frac{\mv{\xi} \cdot \mv{u}}{c_s^2} + \left( \frac{\mv{\xi}\cdot\mv{u}}{c_s^2} \right)^2
  - \frac{\mv{u}^2}{c_s^2}
  + \frac{1}{2} \left( \frac{\mv{\xi}^2}{c_s^2} - D \right) \right]
  + E f_i^{eq}
  \label{eq:equilibrium_total}
  \,.
\f}


## Computing macroscopic variables
The macroscopic variables are calculated by summing over all PPDFs, 
i.e., the density and velocities by calculating

\f{equation}{
  \rho(\mv{x},t) = \sum_i f_i(\mv{x},t),
  \label{eq:lbdensity}
\f}
\f{equation}{
  \rho(\mv{x},t)\mv{u}(\mv{x},t) = \sum_i \xi_i f_i(\mv{x},t).
  \label{eq:lbvelocity}
\f}

For the basic thermal approach, the temperature is computed by
\f{equation}{
  T(\mv{x},t) = \sum_i g_i(\mv{x},t).
  \label{eq:lbtemperature}
\f}

For the approaches considering the internal and total energy, the temperature is derived from
\f{equation}{
  \rho(\mv{x},t)\epsilon = \sum_i h_i(\mv{x},t),
  \label{eq:lbtemperature_inner}
\f}
and

\f{equation}{
  \rho(\mv{x},t)E = \sum_i k_i(\mv{x},t).
  \label{eq:lbtemperature_total}
\f}

The static pressure \f$ p_s \f$ is derived from the ideal gas law by
\f{equation}{
  p_s = \rho RT = \rho c_s^2 = \rho/3 
  \label{eq:lbpressure_static}
\f}


## PPDF directions ## 
The table below shows the PPDF directions of the D3Q27 model.

<table>
<!-- <thead> -->
<tr class="header">
<th><p>0</p></th>
<th><p>1</p></th>
<th><p>2</p></th>
<th><p>3</p></th>
<th><p>4</p></th>
<th><p>5</p></th>
<th><p>6</p></th>
<th><p>7</p></th>
<th><p>8</p></th>
<th><p>9</p></th>
<th><p>10</p></th>
<th><p>11</p></th>
<th><p>12</p></th>
<th><p>13</p></th>
<th><p>14</p></th>
<th><p>15</p></th>
<th><p>16</p></th>
<th><p>17</p></th>
<th><p>18</p></th>
<th><p>19</p></th>
<th><p>20</p></th>
<th><p>21</p></th>
<th><p>22</p></th>
<th><p>23</p></th>
<th><p>24</p></th>
<th><p>25</p></th>
<th><p>26</p></th>
</tr>
<!-- </thead> -->
<!-- <tbody> -->
<tr class="odd">
<td><p>-1<br />
0<br />
0</p></td>
<td><p>1<br />
0<br />
0</p></td>
<td><p>0<br />
-1<br />
0</p></td>
<td><p>0<br />
1<br />
0</p></td>
<td><p>0<br />
0<br />
-1</p></td>
<td><p>0<br />
0<br />
1</p></td>
<td><p>-1<br />
-1<br />
0</p></td>
<td><p>-1<br />
1<br />
0</p></td>
<td><p>1<br />
-1<br />
0</p></td>
<td><p>1<br />
1<br />
0</p></td>
<td><p>-1<br />
0<br />
-1</p></td>
<td><p>-1<br />
0<br />
1</p></td>
<td><p>1<br />
0<br />
-1</p></td>
<td><p>1<br />
0<br />
1</p></td>
<td><p>0<br />
-1<br />
-1</p></td>
<td><p>0<br />
-1<br />
1</p></td>
<td><p>0<br />
1<br />
-1</p></td>
<td><p>0<br />
1<br />
1</p></td>
<td><p>-1<br />
-1<br />
-1</p></td>
<td><p>-1<br />
-1<br />
1</p></td>
<td><p>-1<br />
1<br />
-1</p></td>
<td><p>-1<br />
1<br />
1</p></td>
<td><p>1<br />
-1<br />
-1</p></td>
<td><p>1<br />
-1<br />
1</p></td>
<td><p>1<br />
1<br />
-1</p></td>
<td><p>1<br />
1<br />
1</p></td>
<td><p>0<br />
0<br />
0</p></td>
</tr>
<!-- </tbody> -->
</table>

## Further collision steps

Besides the BGK approach presented above, further collision models are implemented
in the LB solver ([List](@ref collisionStep)). Some of the implemented collision steps
are an combination of the normal BGK collision step with a thermal LB extension. 
Furthermore, the user can choose between incompressible and compressible collision steps,
which are still based on the BGK model. The remaining approaches are more complex
than the BGK model, but might offers advantages for specific applications. For instance,
the MRT model conducts the collision in momentum space. That is, it provides a separate
relaxation frequency for each moment. Hence, it is also referred to as multi-relaxation
model. Its collision equation reads
\f{equation}
  f_i(\mv{x} + \mv{\xi}_i \Delta t, t + \Delta t) = f_i(\mv{x}, t) + \mv{\mv{M}}^{-1}
  \mv{\mv{S}} (\mv{m}^{eq}(\mv{x}, t)-\mv{m}(\mv{x}, t)) \,.
  \label{eq:MRT}
\f}
Due to the use of one relaxation frequency per moment, the MRT provides increased
numerical stability for many simulations. However, the correct adjustment of the
relaxation frequencies might be difficult.

The cascaded LB [Geier2006] and the cumulant [Geier2015] LB are successors of the MRT collision step.
The cascased LB performs the collision using the central moments. The cumulant performs
the collision in cumulant space. These transformation allows the user to simulate high
Reynolds number flows, where the relaxation frequency \f$ \omega \f$ is close to
the limit of \f$ 2 \f$.

In addition to the so far mentioned collision steps, nearly all models are equipped with
a Smagorinsky model to perform an LES. To activte the Smagorinsky model the property
`solverMethod` must be adapted.

## Non-dimensional Variables {#LB_nonDim}

For the calculation non-dimensional variables and parameters are used.
The reference values are the density @f$\rho_0@f$, the mean thermodynamic
velocity @f$\xi_0@f$, the grid distance @f$\Delta x@f$ and the temperature
@f$T_0@f$. In the code implementation the following non-dimensional
variables are used:
<table>
<tr><th>Variable</th>                     <th>Equation</th></tr>
<tr><td>Velocity                     </td>  <td>@f$u^*=\frac{u}{\xi_0}@f$</td></tr>
<tr><td>Isothermal speed of sound    </td>  <td>@f$c_s^*=\frac{c_s}{\xi_0}=\sqrt{\frac{1}{3}}@f$</td></tr>
<tr><td>Molecular velocity           </td>  <td>@f$\xi_0^*=\frac{\xi_0}{\xi_0}=1@f$</td></tr>
<tr><td>Grid distance                </td>  <td>@f$\Delta x^*=\frac{\Delta x}{\Delta x}=1@f$</td></tr>
<tr><td>Characteristic Length        </td>  <td>@f$L^*=\frac{L}{\delta x}@f$</td></tr>
<tr><td>Time step                    </td>  <td>@f$\Delta t^*=\frac{\Delta t \cdot \xi_0}{\Delta x} = 1@f$</td></tr>
<tr><td>Density                      </td>  <td>@f$\rho^*=\frac{\rho}{\rho_0}@f$</td></tr>
<tr><td>Pressure                     </td>  <td>@f$p^*=\frac{p}{\rho_0 \cdot \xi_0^2}@f$</td></tr>
<tr><td>Collision Operator           </td>  <td>@f$\Omega^* = \omega \cdot \Delta t = \frac{2}{6 \cdot \nu^* + 1}@f$</td></tr>
<tr><td>Temperature                  </td>  <td>@f$T^*=\frac{T}{T_0}@f$</td></tr>
<tr><td>Thermal diffusivity          </td>  <td>@f$\kappa^*=\frac{\nu^*}{Pr}@f$</td></tr>
<tr><td>Collision Operator (thermal) </td>  <td>@f$\Omega_g^* = \omega_g \cdot \Delta t = \frac{2}{6 \cdot \frac{\nu^*}{Pr} + 1}@f$</td></tr>
</table>

The non-dimensional kinematic viscosity @f$\nu^*@f$ can be calculated using
the Reynolds similarity (@f$Re = Re^*@f$) and the Mach-Number
@f$Ma^* = \frac{u^*}{c_s^*} = \sqrt{3} \cdot u^*@f$ (which is given in the
properties_run file)

@f{equation}{Re = \frac{L \cdot u}{\nu} = \frac{L^* \cdot u^*}{\nu^*} \Leftrightarrow \nu^* = \frac{u^* \cdot L^*}{Re} = \frac{L}{\sqrt{3} \cdot \Delta x} \cdot \frac{Ma^*}{Re}@f}

Furthermore, it can be shown, that the Mach-Number @f$Ma^*@f$ only affects
the temporal accuracy:

@f{equation}{Ma^* = \frac{u^*}{c_s^*} = \sqrt{3} \cdot u^* = \sqrt{3} \cdot \frac{u}{\xi_0} \Leftrightarrow \xi_0 = \frac{\sqrt{3} \cdot u}{Ma^*} \Leftrightarrow \Delta t = Ma^* \cdot \frac{\Delta x}{\sqrt{3} \cdot u}@f}

By determining the time step @f$\Delta t@f$ and afterwards the mean
thermodynamic velocity @f$\xi_0 = \frac{\Delta x}{\Delta t}@f$ the velocity
@f$ u @f$ can be calculated. The density @f$\rho@f$ and the temperature @f$ T @f$ can easily
be calculated from the non-dimensional variable and the corresponding
reference value. 


## Boundary Conditions {#boundary_conditions}
In the following, the most important boundary conditions are explained. A list of all
boundary conditions can be found here: @subpage LBBCList

### Inflow boundary conditions

`BC 1000` `BC 1060 (Thermal)`

Velocity is prescribed, the boundary normal is used as velocity vector,
the pressure is extrapolated in opposite direction (von Neumann
condition). Define the properties

` lbControlInflow = 1`

` inflowSegmentIds = 1,2,... ;`

to prescribe a velocity field in a parabolic manner.

`BC 1002` `BC 1022 (Thermal)`

Velocity is prescribes. Use this for a periodic channel, main direction
is X, periodic in Z-direction. channel center (Y-component) is placed at
y=0. Profile is given by ($h$ is half of channel height):

\f$u(y)=u_{max}\left(h^2-y^2\right)=\frac{3}{2}\frac{\bar{u}}{h^2}\left(h^2-y^2\right)=\frac{3}{2}\frac{\frac{Ma}{\sqrt{3}}}{h^2}\left(h^2-y^2\right)\f$

`BC 1061 (only thermal)`

Prescribes temperature and velocity with a Blasius profile for
boundary-layer flows. 
@warning Works only in serial at the moment!

`BC 1080`

Generates a set of boundary conditions for in- and outflows of internal
flow configurations such as a pipe. The temporal variation of the volume
flow (or the main stream velocity) is reconstructed using Fourier
coefficients determined by the FFTW library. The non-dimensional
frequency \f$\alpha\f$ (m_alpha) is defined as

@f{equation}{\alpha = \frac{\mathrm{Wo}}{\sqrt{2\pi}} = \frac{\mathtt{m\_diameter}}{\sqrt{\mathtt{m\_nu}\times N}}@f}

where \f$N\f$ is the number of time-steps per one cycle and the Womersley
number is \f$\mathrm{Wo}=L\sqrt{\omega/\nu}\f$ at an angular frequency
\f$\omega\f$. For a given Reynolds number $\mathrm{Re}=UL/\nu$ and a Mach
number \f$\mathrm{Ma}=U/\sqrt{RT}\f$ the number of time steps is

@f{equation}{N = \frac{2\pi \mathrm{Re}}{\mathrm{Wo}^2}\times \mathtt{m\_diameter} \times \frac{\sqrt{3}}{\mathrm{Ma}}@f}

In practice a validation can use the flow configuration in literature
(Varghese et al., J. Fluids Mech., Vol. 582, 2007). In September 2019
the simulation using LB (120 and 200 million cells) shows a good
agreement with the DNS results of stenotic flows at a Reynolds number
600.

### Wall conditions (no slip) {#wall_conditions_no_slip}

`BC 2000`

interpolated bounce back with BFL-rule of Bouzidi M, Firdaouss M,
Lallemand P. Momentum transfer of a Boltzmann-lattice fluid with
boundaries. Physics of Fluids.
2001;13(11):3452-3459.

`BC 2001`

simple bounce back, no interpolation, wall is located at the cell
surface

`BC 2002`

simple bounce back, no interpolation, wall is located at the cell center

`BC 2020`

Same as BC 2000, but also sets the temperature to the value specified by
the property
```
initTemperatureKelvin.
```

### Slip conditions {#slip_conditions}

For a slip condition specular reflection of the PPDFs is performed at
the surface.

`BC 302d`

Replace d by the direction for which slip should be implemented

      +-x   +-y   +-z
  --- ----- ----- -----
  d   0,1   2,3   4,5

### Zero gradient and Neumann type conditions

`BC 3000`

Zero gradient in direction of boundary normal

`BC 301d`

Zero gradient in axes-directions (d: see below)

`<b>`{=html}Definition of d:`</b>`{=html} the extrapolation direction is
always pointing out of the domain

      -x   +x   -y   +y   -z   +z
  --- ---- ---- ---- ---- ---- ----
  d   0    1    2    3    4    5


### Outflow conditions

`BC 4000` `BC 4110 (thermal)`

Pressure is prescribed, the velocity is extrapolated normal to the
surface

`BC 4071` `BC 4081 (thermal)`

Pressure is prescribed, the velocity is extrapolated normal to the
surface. The prescribed pressure is read from the property
`<i>`{=html}`<b>`{=html}rho1`</b>`{=html}`</i>`{=html}.

### Outflow condition for variable pressure

`BC 4070` `BC 4080`

Pressure is prescribed according to appearing volume flux and pressure
distribution from last iteration (use `BC 4080` for thermal simulations)

Using the equation of Saint Vernant and Wanzel, one can obtain the
pressure by using

@f$p_1 = p_0\left(1 - \frac{\gamma - 1}{2\gamma}\cdot\frac{\rho_0}{p_0}\cdot v_0^2\right)^{\frac{\gamma}{\gamma -1}}@f$

A non-dimensionalization with

@f$\rho^\ast = \frac{\rho}{\rho_0},\quad c_{s}^{\ast} = \frac{c_s}{\xi_0},\quad v^{\ast} = \frac{v}{\xi_0}@f$

and using $p=\rho c_s^2$, $c_s=\frac{1}{\sqrt{3}}\xi_0$ and extending
with $\frac{\rho^{\ast^2}}{\rho^{\ast^2}}$ leads to

@f$\rho_{i+1}=\left(1 - \frac{\gamma - 1}{2\gamma}\cdot \frac{3}{\rho^{\ast^2}_{i}}\left(\rho^{\ast} v^\ast\right)^2\right)^{\frac{\gamma}{\gamma -1}}@f$

So, first all the velocities are extrapolatet to the boundary cells and the new density ratio is calculated.

`BC 4072` `BC 4082`

Pressure is prescribed according to a Reynolds number that should be
reached. Define the density via the property
```
\f$\rho_1\f$
```
\f$\delta\rho = 1.0 - \rho_1,\f$ which is used for the maximal step size
to in- or decrease the density ratio. In each timestep the local
Reynolds number is calculated and the desity ratio to be prescribed is
adjusted accordingly.

`BC 4073`

BC 4073 works similar as BC 4072, except that the Reynolds number is measured at a
specified location which can be defined in the geometry file by
```
localReCut = P1, P2, P3, N1, N2, N3

localReCutInterval = interval (defines the interval to adapt the density - it requires some time to propagate the information)

localReCutReportInterval = report interval (when to write to the file BCResidual.dat)

localReCutDistance = distance (the radius in the plane for the definition of the cells to consider for the local Re calculation)

localReCutAdpPerc = percentage (0.0 - 1.0) (when to start adaption based on the difference of the actual value to the desired value of the Re number)
```
where P1, P2, P3, N1, N2, and N3 define a plane point and normal.

If you are running a restart define:
```
localReCut_rho1 = density for initialization

localReCut_lRho = local density
 
localReCut_rhoLast = density from last iteration

localReCut_deltaRho = delta for the density

localReCut_ReLast = last adaption Re number
```

@note Note that these values can be obtained from your BCResidual.dat file.

### Periodic Boundary Conditions {#periodic_boundary_conditions}

To generate a perodic link the corresponding stl files are declared as
periodic segment in the geometry file. No bc is needed.

For example:
```
geometry.toml

periodicSegmentIds = [0,1]

filename.0 = "stl/face_x_neg.stl"

filename.1 = "stl/face_x_pos.stl"

bc.0 = 0

bc.1 = 0

properties_run.toml

periodicCartesianDir = [1, 0, 0]
```

## Flow Field Initialization {#flow_field_initialization}

The property `initMethods` determines the LBM initalization method.
The LB solver provides several initialization conditions. Some of these are quite specific. Others are more general.
In the following, only the general ICs are presented. A full list of ICs can be found here: @subpage LBICList
```
LB_FROM_ZERO_INIT (sets velocity to zero and density and temperature to unity)

LB_LAMINAR_INIT_PX (sets velocity in main flow direction (here x-direction))

LB_LAMINAR_PIPE_INIT (sets parabolic velocity profil in main flow direction)

LB_LAMINAR_CHANNEL_INIT (sets parabolic velocity profile for channel like geometries in main flow direction)
```

## References

* H. Chen, S. Chen, and W. H. Matthaeus, Recovery of the Navier-Stokes equations
using a lattice-gas Boltzmann method, Physical Review A 45.8, 1992, R5339–
R5342, [10.1103/PhysRevA.45.R5339][Chen1992].

[Chen1992]: https://doi.org/10.1103/PhysRevA.45.R5339

* Y. H. Qian, D. D’Humières, and P. Lallemand, Lattice BGK Models for Navier-
Stokes Equation, Europhysics Letters 17.6, 1992, p. 479, [10.1209/0295-5075/17/6/001][Qian1992].

[Qian1992]: https://doi.org/10.1209/0295-5075/17/6/001

* X. Shan, Simulation of Rayleigh-Benard Convection Using a Lattice Boltzmann
Method, Physical Review E 55.3, 1997, pp. 2780–2788, [10.1103/PhysRevE.55.2780][Shan1997]. 

[Shan1997]: https://doi.org/10.1103/PhysRevE.55.2780

* X. He, S. Chen, and G. D. Doolen, A Novel Thermal Model for the Lattice Boltzmann
Method in Incompressible Limit, Journal of Computational Physics
146.1, 1998, pp. 282–300, [10.1006/jcph.1998.6057][He1998].

[He1998]: https://doi.org/10.1006/jcph.1998.6057

* Z. Guo, C. Zheng, B. Shi, and T. S. Zhao, Thermal lattice Boltzmann equation
for low Mach number flows: Decoupling model, Physical Review E 75.3, 2007,
p. 036704, [10.1103/PhysRevE.75.036704][Guo2007].

[Guo2007]: https://doi.org/10.1103/PhysRevE.75.036704

* M. Geier, A. Greiner, and J. G. Korvink, Cascaded digital lattice Boltzmann automata
for high Reynolds number flow, Physical Review E 73, 2006, p. 066705,
[10.1103/PhysRevE.73.066705][Geier2006]

[Geier2006]: https://doi.org/10.1103/PhysRevE.73.066705

* M. Geier, M. Schoenherr, A. Pasquali, M. Kraftczyk, The cumulant lattice Boltzmann
equation in three dimensions: Theory and validation, Computers \& Mathematics with
Applications 70.4, 2015, 507-547, [10.1016/j.camwa.2015.05.001][Geier2015]

[Geier2015]: https://doi.org/10.1016/j.camwa.2015.05.001

-   Bhatnagar PL, Gross EP, Krook M. A Model for Collision Processes in
    Gases. I. Small Amplitude Processes in Charged and Neutral
    One-Component Systems. Phys. Rev. 1954;94(3):511-525.
-   Chen H, Chen S, Matthaeus WH. Recovery of the Navier-Stokes
    equations using a lattice-gas Boltzmann method. Phys. Rev. A.
    1992;45(8):R5339\--R5342.
-   Chopard B, Droz M. Cellular Automata Modeling of Physical Systems.
    Cambridge: Cambridge University Press; 1998.\
    Available at: <http://ebooks.cambridge.org/ref/id/CBO9780511549755>.
-   Dubois F, Lallemand P. Towards higher order lattice Boltzmann
    schemes. Journal of Statistical Mechanics: Theory and Experiment.
    2009(06):P06006.\
    Available at:
    <http://stacks.iop.org/1742-5468/2009/i=06/a=P06006?key=crossref.e63ddccc4fcebce54fb855d8e1f53e07>.
-   Eitel G, Freitas R, Lintermann A, Meinke M, Schröder W. Numerical
    Simulation of Nasal Cavity Flow Based on a Lattice-Boltzmann Method.
    In: Dillmann A, Heller G, Klaas M, et al., eds. New Results in
    Numerical and Experimental Fluid Mechanics VII.Vol 112. Springer
    Berlin / Heidelberg; 2010:513-520.
-   Eitel G, Soodt T, Schröder W. Investigation of pulsatile flow in the
    upper human airways. International Journal of Design & Nature and
    Ecodynamics. 2010;5(4).
-   Hartmann D, Meinke M, Schröder W. A General Formulation of Boundary
    Conditions on Cartesian Cut-Cells for Compressible Viscous Flow. In:
    19th AIAA CFD conference. San Antonio, TX: American Institute of
    Aeronautics and Astronautics, 1801 Alexander Bell Dr., Suite 500
    Reston VA 20191-4344 USA,; 2009:1-13.\
    Available at:
    <http://pdf.aiaa.org/preview/2010/CDReadyMCFD09_2126/PV2009_3878.pdf>.
-   Hänel D. Molekulare Gasdynamik. Springer-Verlag; 2004.
-   Nourgaliev R. The lattice Boltzmann equation method: theoretical
    interpretation, numerics and implications. International Journal of
    Multiphase Flow. 2003;29(1):117-169.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0301932202001088>.
-   Qian YH, D'Humières D, Lallemand P. Lattice BGK Models for
    Navier-Stokes Equation. Europhysics Letters (EPL).
    1992;17(6):479-484.\
    Available at:
    <http://stacks.iop.org/0295-5075/17/i=6/a=001?key=crossref.9c570b56da56c5c15421217daa9bb1dd>.
-   Shu C, Niu XD, Chew YT, Cai QD. A fractional step lattice Boltzmann
    method for simulating high Reynolds number flows. Mathematics and
    Computers in Simulation. 2006;72(2-6):201-205.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0378475406001480>.
-   Succi S. The Lattice Boltzmann equation: theory and applications.
    Rome: Oxford University Press; 2001.
-   Wagner AJ. A Practical Introduction to the Lattice Boltzmann Method.
    North. 2008;(March).
-   Wolf-Gladrow DA. Lattice-Gas Cellular Automata and Lattice Boltzmann
    Models - An Introduction. Springer Verlag; 2000:308. Available at:\
    <http://books.google.com/books?hl=en&lr=&id=cHcpWgxAUu8C&oi=fnd&pg=PA15&dq=Lattice-Gas+Cellular+Automata+and+Lattice+Boltzmann+Models+-+An+Introduction&ots=I7s4HstRrd&sig=n1vEIVSOofn4XDik1p_R-r9zsK8>.
-   Xiu-Ying K, Yu-Pin J, Da-He L, Yong-Juan J. Three-dimensional
    lattice Boltzmann method for simulating blood flow in aortic arch.
    Chinese Physics B. 2008;17(3):1041-1049.\
    Available at:
    <http://stacks.iop.org/1674-1056/17/i=3/a=049?key=crossref.75ec865996830dffc9c7ceff0fe91e52>.

### Boundary conditions {#boundary_conditions_1}

-   Bao J, Yuan P, Schaefer L. A mass conserving boundary condition for
    the lattice Boltzmann equation method. Journal of Computational
    Physics. 2008;227(18):8472-8487.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0021999108003288>.
-   Bouzidi M, Firdaouss M, Lallemand P. Momentum transfer of a
    Boltzmann-lattice fluid with boundaries. Physics of Fluids.
    2001;13(11):3452-3459.\
    Available at:
    <http://link.aip.org/link/PHFLE6/v13/i11/p3452/s1&Agg=doi>.
-   Caiazzo A, Junk M. Boundary forces in lattice Boltzmann: Analysis of
    Momentum Exchange algorithm. Computers & Mathematics with
    Applications. 2008;55(7):1415-1423.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0898122107006281>.
-   Caiazzo A, Maddu S. Lattice Boltzmann boundary conditions via
    singular forces: Irregular expansion analysis and numerical
    investigations. Computers & Mathematics with Applications.
    2009;58(5):930-939.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0898122109000984>.
-   Chen S, Martínez D, Mei R. On boundary conditions in lattice
    Boltzmann methods. Physics of Fluids. 1996;8(9):2527.\
    Available at:
    <http://link.aip.org/link/PHFLE6/v8/i9/p2527/s1&Agg=doi>.
-   Geveler M, Ribbrock D, Göddeke D, Turek S. Lattice-Boltzmann
    Simulation of the Shallow-Water Equations with Fluid-Structure
    Interaction on Multi-and Manycore Processors. In: Facing the
    Multicore Challenge. Heidelberg: Springer LNCS; 2010.\
    Available at:
    <http://www.mathematik.tu-dortmund.de/~goeddeke/pubs/pdf/Geveler_2010_LBS.pdf>.
-   Guo Z, Zheng C, Shi B. An extrapolation method for boundary
    conditions in lattice Boltzmann method. Physics of Fluids.
    2002;14(6):2007.\
    Available at:
    <http://link.aip.org/link/PHFLE6/v14/i6/p2007/s1&Agg=doi>.
-   Junk M, Yang Z. Outflow boundary conditions for the lattice
    Boltzmann method. Progress in Computational Fluid Dynamics, an
    International Journal. 2008;8(1):38-48.\
    Available at:
    <http://inderscience.metapress.com/index/y27727h837111461.pdf>.
-   Junk M, Yang Z. Asymptotic Analysis of Lattice Boltzmann Boundary
    Conditions. Journal of Statistical Physics. 2005;121(1-2):3-35.\
    Available at:
    <http://www.springerlink.com/index/10.1007/s10955-005-8321-2>.
-   Maier RS, Bernard RS, Grunau DW. Boundary conditions for the lattice
    Boltzmann method. Physics of Fluids. 1996;8(7):1788.\
    Available at:
    <http://link.aip.org/link/PHFLE6/v8/i7/p1788/s1&Agg=doi>.
-   Mei R, Luo L-S, Shyy W. An Accurate Curved Boundary Treatment in the
    Lattice Boltzmann Method. Journal of Computational Physics.
    1999;155(2):307-330.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0021999199963349>.
-   Mei R, Shyy W, Yu D. Lattice Boltzmann Method for 3-D Flows with
    Curved Boundary. Journal of Computational Physics.
    2000;161(2):680-699.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0021999100965227>.
-   Ziegler DP. Boundary conditions for lattice Boltzmann simulations.
    Journal of Statistical Physics. 1993;71(5-6):1171-1177.\
    Available at:
    <http://www.springerlink.com/index/10.1007/BF01049965>.
-   Zou Q, He X. On pressure and velocity boundary conditions for the
    lattice Boltzmann BGK model. Physics of Fluids. 1997;9(6):1591.\
    Available at:
    <http://link.aip.org/link/PHFLE6/v9/i6/p1591/s1&Agg=doi>.

### Heat transfer {#heat_transfer}

-   Chew YT, Niu XD, Shu C. Three-dimensional lattice Boltzmann BGK
    model and its application to flows with heat transfer in a
    rectangular microchannel. International Journal for Numerical
    Methods in Fluids. 2006;50(11):1321-1334.\
    Available at: <http://doi.wiley.com/10.1002/fld.1143>.
-   Guo Z, Shi B, Zheng C. A coupled lattice BGK model for the
    Boussinesq equations. International Journal for Numerical Methods in
    Fluids. 2002;39(4):325-342.\
    Available at: <http://doi.wiley.com/10.1002/fld.337>.
-   He X, Chen S, Doolen GD. A Novel Thermal Model for the Lattice
    Boltzmann Method in Incompressible Limit. Journal of Computational
    Physics. 1998;146(1):282--300.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0021999198960570>.
-   Meng F, Wang M, Li Z. Lattice Boltzmann simulations of oscillating
    flow and heat transfer under linearization assumption. In:
    Proceeding of the First International Conference on Enhancement and
    Promotion of Computational Methods in Engineering Science and
    Mechanics, China. Beijing, China; 2006:474--478.\
    Available at: <http://cnls.lanl.gov/~mwang/pdf/2006JLU.pdf>.
-   Olander J. COMPARISON OF THE HYBRID AND THERMAL LATTICE-BOLTZMANN
    METHODS. 2009;(August).\
    Available at: <http://smartech.gatech.edu/handle/1853/31826>.
-   Peng Y, Shu C, Chew YT. A 3D incompressible thermal lattice
    Boltzmann model and its application to simulate natural convection
    in a cubic cavity. Journal of Computational Physics.
    2004;193(1):260--274.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0021999103004303>.
-   Shi Y, Zhao T, Guo Z. Thermal lattice Bhatnagar-Gross-Krook model
    for flows with viscous heat dissipation in the incompressible limit.
    Physical Review E. 2004;70(6):1-10.\
    Available at: <http://link.aps.org/doi/10.1103/PhysRevE.70.066310>.
-   Vahala G, Pavlo P, Vahala L, Martys NS. Thermal Lattice-Boltzmann
    Model (TLBM) for Compressible Flows. International Journal of Modern
    Physics C-Physics and Computer. 1998;9(8):1247--1262.\
    Available at:
    <http://www.worldscinet.com/abstract?id=pii:S0129183198001126>.
-   Weiss JP. Numerical analysis of lattice Boltzmann methods for the
    heat equation on a bounded interval. Institut für Angewandte und
    Numerische Mathematik (Inst. f. Angew. u. Num. Math.),
    Universitätsverlag Karlsruhe, Karlsruhe; 2006.\
    Available at: <http://nbn-resolving.de/urn:nbn:de:0072-53048>.
-   Zheng L, Shi BC, Cai Z. Lattice Boltzmann Method for Simulating the
    Temperature Jump and Velocity Slip in Microchannels. Communications
    in Computational Physics. v2. 2007;1125(6):1125-1138.\
    Available at:
    <http://cscs.hust.edu.cn/uploadfile/download/guozl/ZhengL_CCP2007_>
    Lattice Boltzmann Method for Simulating the Temperature Jump
    andVelocity Slip inMicrochannels.pdf.

### LES

-   Menon S, Soo J. Direct and large-eddy simulations of
    three-dimensional jets using the lattice Boltzmann method. Aerospace
    Engineering. 2002;(September).\
    Available at: <http://smartech.gatech.edu/handle/1853/12013>.
-   YU H, GIRIMAJI S, LUO L. DNS and LES of decaying isotropic
    turbulence with and without frame rotation using lattice Boltzmann
    method. Journal of Computational Physics. 2005;209(2):599-616.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0021999105001907>.

### Refinement

-   Dupuis A, Chopard B. Theory and applications of an alternative
    lattice Boltzmann grid refinement algorithm. Physical Review E.
    2003;67(6):1-7.\
    Available at: <http://link.aps.org/doi/10.1103/PhysRevE.67.066707>.
-   Filippova O. Grid Refinement for Lattice-BGK Models. Journal of
    Computational Physics. 1998;147(1):219-228.\
    Available at:
    <http://linkinghub.elsevier.com/retrieve/pii/S0021999198960892>.
-   Filippova O, Hänel D. Boundary-Fitting and Local Grid Refinement for
    Lattice-BGK Models. International Journal of Modern Physics C.
    1998;9(8):1271-1279.


