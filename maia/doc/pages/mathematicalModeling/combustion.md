# Combustion Equations# {#mmCombustion}

## Detailed chemistry formulation

#### Disclaimer: The equations given here a dimensional!

### Selection of conservative variables set and general equations

A set of primitive and conservative variables is defined for two-dimensional flows. Note that the extension to three-dimensional flows is straightforward. The set of conservative variables, namely \f$\mathrm{G}_c\f$, is the non-reacting Navier-Stokes equations conservative variables set expanded by a total of \f$N\f$ additional species mass fraction variables \f$Y_k\f$, required to define the chemical composition of the gas mixture, where the subscript \f$k\f$ denotes the associated species. The set of primitive variables, \f$G_p\f$, is derived from \f$G_c\f$ and the equation of state. In this way, the conservative variables and primitive variables are chosen as 
$$
\mathrm{G}_c=\left[\rho, \rho \mathbf{U}, \rho E, \rho Y_1, \ldots, \rho Y_N\right]
$$

$$
\mathrm{G}_p=\left[\rho, \mathbf{U}, p, Y_1, \ldots, Y_N\right]
$$

, where \f$\mathbf{U}\f$ is the velocity vector, \f$\rho\f$ defines the density of the fluid, \f$E\f$ is the total energy and $p$ is the pressure. A total of \f$(N+4)\f$ variables are therefore necessary to completely describe the problem and to solve the corresponding equation system for two-dimensional flows.  

The introduced mass fraction \f$Y_k\f$ is defined as

\f{align}{
Y_k=\frac{m_k}{m},
\f}

where \f$m_k\f$ is the mass of species \f$k\f$ in a given volume \f$V\f$, and $m$ is the total mass in the same given volume. The mass fraction is related to the molar fraction \f$X_k\f$, by the following equation

\f{align}{
X_k=\frac{\bar{W}}{W_k} Y_k,
\f}

where \f$\bar{W}\f$ is the mean molar weight of the mixture and \f$W_k\f$ is the species molar weight. The mean molar weight can be computed through

\f{align}{
\frac{1}{\bar{W}}=\sum_{k=1}^N \frac{Y_k}{W_k} \quad \text { or } \quad \bar{W}=\sum_{k=1}^N X_k W_k
\f}

For the energy equation the total non-chemical energy was chosen as the dependent variable. The total energy is related to the temperature field and velocity field through the expression

\f{align}{
E=e_s+\frac{1}{2} \mathbf{U U}=\int_{T_0}^T c_v(T) d T-\frac{R T_0}{\bar{W}}+\frac{1}{2} \mathbf{U U},
\f}

where \f$E\f$ stands for the total non-chemical energy, $e_s$ is the sensible energy, \f$T_0\f$ is the reference temperature and \f$c_v\f$ is the specific heat of the mixture at constant volume. The term \f$R T_0 / \bar{W}\f$ is the sensible energy level at the reference temperature \f$T_0\f$, which is unequal zero. The sensible energy is related to the other forms of energy through the expressions given in the following table:

| Form | Energy | Enthalpy |
| :--- | :--- | :--- |
| Sensible | \f$e_s=\int_{T_0}^T c_v d T-R T_0 / \bar{W}\f$ | \f$h_s=\int_{T_0}^T c_p d T\f$ |
| Sensible \f$+\f$ Chemical | \f$e=e_s+\sum_{k=1}^N \Delta h_{f, k}^o Y_k\f$ | \f$h=h_s+\sum_{k=1}^N \Delta h_{f, k}^o Y_k\f$ |
| Total Chemical | \f$e_t=e+\frac{1}{2} \mathbf{U U}\f$ | \f$h_t=h+\frac{1}{2} \mathbf{U U}\f$ |
| Total non Chemical | \f$E=e_s+\frac{1}{2} \mathbf{U U}\f$ | \f$H=h_s+\frac{1}{2} \mathbf{U U}\f$ |

A relation can be derived from the first law of thermodynamics to relate both heat capacities \f$c_p\f$ and \f$c_v\f$ for an ideal gas mixture. Similarly, an expression can be derived for the pure species heat capacities \f$c_{p, k}\f$ and \f$c_{v, k}\f$, which are required to compute the pure species sensible enthalpies. Both relations can be written as

$$
c_p-c_v=\frac{R}{\bar{W}} \quad \text { and } \quad c_{p, k}-c_{v, k}=\frac{R}{W_k},
$$

where \f$R\f$ is the universal gas constant and \f$c_p\f$ is the specific heat of the mixture at constant pressure. To completely close the set of equations, the pressure is computed from the temperature and density fields through an equation of state (ideal gas law), expressed as

$$
p=\frac{\rho T R}{\bar{W}} \text {. }
$$


### Heat capacity computation
The heat capacity is a species thermodynamic property and depends both on the temperature and on the pressure. However, the pressure dependence is usually small compared to the temperature dependence and is usually ignored. No exact, analytical expressions exist to compute the heat capacities as a function of the temperature, however, approximate expressions were developed by the NASA in the 90s. The heat capacity is approximated by a fourth-order polynomial, with coefficients \f$a_i, i \in[0,1,2,3,4]\f$, given as thermodynamic data to the chemistry solver. The polynomial has the form

$$
W_k \frac{c_{p, k}}{R}=a_{k, 0}+a_{k, 1} T+a_{k, 2} T^2+a_{k, 3} T^3+a_{k, 4} T^4 .
$$

A similar expression follows for \f$c_{v, k}\f$
$$
W_k \frac{c_{v, k}}{R}=\left(a_{k, 0}-1\right)+a_{k, 1} T+a_{k, 2} T^2+a_{k, 3} T^3+a_{k, 4} T^4 .
$$

Due to limitations in the fitting procedure for the coefficients \f$a_i\f$, two different sets of coefficients exist for the low temperature and high temperature region, \f$a_{i, l o w}\f$ and \f$a_{i, \text { high }}\f$, respectively. In addition, the mixture specific heat capacity is computed as a mass weighted summation of the individual species specific heat capacities as
$$
c_p=\sum_{k=1}^N Y_k c_{p, k} .
$$

The analytical expression for the total energy after integrating over the low temperature region reads
$$
E=\sum_{k=1}^N\left(\frac{Y_k R}{W_k}\left(\sum_{i=0}^4 \frac{a_{k, i, l o w}}{i+1} T^{* i+1}-T^*\right)_{T_0}^T\right)-\frac{R T_0}{\bar{W}}+\frac{1}{2} \mathbf{U U}
$$

### Transport equations

#### Species conservation equations
The mass conservation equation for species \f$k\f$ can be written as

\f{align}{
\frac{\partial \rho Y_k}{\partial t}+\frac{\partial}{\partial x_i}\left(\rho\left(u_i+V_{k, i}\right) Y_k\right)=\dot{\omega}_k,
\f}

where \f$\dot{\omega}_k\f$ and \f$V_k\f$ are the reaction rate and the diffusion velocity of species $k$ respectively. The term \f$\rho Y_k V_k\f$ is called species mass flux \f$j_k\f$. The species diffusive velocity \f$V_k\f$ can be computed through different methods.  

The following constraints apply for the species conservation equations

\f{align}{
\sum_{k=1}^N j_k=0 \quad \text { and } \quad \sum_{k=1}^N Y_k=1 \quad \text { and } \quad \sum_{k=1}^N \dot{\omega}_k=0 .
\f}

#### Mass conservation equation
The total mass conservation equation remains unchanged compared to non-reacting flows. This is expected since mass is neither created nor destroyed during the chemical reaction process. The mass conservation equation reads

\f{align}{
\frac{\partial \rho}{\partial t}+\frac{\partial\left(\rho u_i\right)}{\partial x_i}=0 .
\f}

Since the mass conservation equation can be algebraically derived from the species transport equations, it is not an independent equation. Consequently, the set of differential equations is overdetermined and only \f$(N-1)\f$ additional species conservation equations need to be solved together with the mass conservation equation.  

Alternatively, the \f$N\f$ species conservation equations can be solved. The system's overdetermination poses no additional challenge when dealing with exact expressions for the species diffusion velocity since the mass conservation is ensured. However, this is no longer the case when first-order approximations are used to approximate the species diffusion velocities, and some additional measures are necessary to ensure global mass conservation.

#### Momentum conservation equation
The momentum conservation equation remains unchanged compared to non-reacting flows:

\f{align}{
\frac{\partial}{\partial t} \rho u_j+\frac{\partial}{\partial x_i} \rho u_i u_j=-\frac{\partial p}{\partial x_j}+\frac{\tau_{i j}}{\partial x_i}+\rho \sum_{k=1}^N Y_k f_{k, j}
\f}

where $f_{k, j}$ are the body forces in direction $i$. 

#### Energy conservation equation

There are multiple forms for the energy conservation equation. The total nonchemical energy \f$E\f$ is chosen as the dependent variable for the formulation found in @maia. The energy conservation equation is therefore written as

\f{align}{
\frac{\partial \rho E}{\partial t}+\frac{\partial \rho u_i E}{\partial x_i}=\dot{\omega}_T-\frac{\partial q_i}{\partial x_i}+\frac{\partial}{\partial x_j}\left(\sigma_{i j} u_i\right)+\dot{Q}+\rho \sum_{k=1}^N Y_k f_{k, j}\left(u_i+V_{k, i}\right).
\f}

The heat release term \f$\dot{\omega}_T\f$ is computed as

\f{align}{
\dot{\omega}_T=-\sum_{k=1}^N \Delta h_{f, k}^o \dot{\omega}_k,
\f}

where \f$\Delta h_{f, k}^o\f$ is the standard formation enthalpy of species \f$k\f$ and \f$\dot{\omega}_k\f$ is the species reaction rate. The standard formation enthalpy is defined as the standard-state enthalpy \f$h^o\f$ minus the standard-state enthalpy of the elements forming the molecular compound. The superscript indicates standard-state conditions, which are taken to be a pressure of 1 bar. The standard temperature can be chosen freely. However, a single reference temperature has to be defined to ensure the consistency of the results. The reference temperature is \f$T_0=298 \mathrm{~K}\f$. From its definition, it follows that for a pure species molecule at the standard pressure and temperature \f$\Delta h_f^o=0\f$. The energy flux \f$q_i\f$ is defined as

\f{align}{
q_i=-\lambda \frac{\partial T}{\partial x_i}+\rho \sum_{k=1}^N h_{s, k} Y_k V_{k, i}-\sum_{k=1}^N \frac{R T}{W_k X_k} D_k^T d_{k, i}
\f}

where \f$\lambda\f$ is the heat transfer coefficient, \f$h_{s, k}\f$ is the species sensible enthalpy, \f$D_k^T\f$ is the thermal diffusivity of species \f$k\f$ and \f$d_k\f$ is the species diffusion driving force. The species' sensible enthalpy is computed as

\f{align}{
h_{s, k}=\int_{T_0}^T c_{p, k}(T) d T .
\f}

Special attention is given to the additional transport terms that appear in the conservation equations. These transport
terms are all diffusive in nature, as they involve the second derivative of a field variable. This, in turn, implies that the convective flux remains unchanged, and only the diffusive flux is affected by the additional terms introduced by the species transport.

### Species mass flux
Before proceeding with the discussion of the different models for the computation of the transport coefficients, the species mass diffusion velocity \f$V_k\f$ is reintroduced, defined as

\f{align}{
\mathbf{V}_{\mathbf{k}}=\mathbf{V}-\tilde{\mathbf{V}}_{\mathbf{k}}
\f}

where \f$\mathbf{V}\f$ is the mass average velocity and \f$\tilde{\mathbf{V}}_{\mathbf{k}}\f$ is the average velocity of species \f$k\f$ relative to the fixed frame of reference. The mass average velocity is defined as

\f{align}{
\mathbf{V}=\frac{1}{\rho} \sum_{k=1}^N \rho_k \tilde{\mathbf{V}}_{\mathbf{k}}=\sum_{k=1}^N Y_k \tilde{\mathbf{V}}_{\mathbf{k}} .
\f}

Once the mass diffusion velocity \f$V_k\f$ is computed, the species mass flux \f$j_k\f$ can be computed as

\f{align}{
j_k=\rho_k\left(\tilde{\mathbf{V}}_{\mathbf{k}}-\mathbf{V}\right)=\rho_k V_k=\rho Y_k V_k .
\f}

Similar expressions can be derived for molar quantities. The molar flux of species \f$k\f$ related to the molar averaged velocity, defined as \f$V^*\f$, is given by

\f{align}{
J_k^*=\left[X_k\right]\left(\tilde{\mathbf{V}}_{\mathbf{k}}-\mathbf{V}^*\right)=\left[X_k\right] \mathbf{V}_{\mathbf{k}}^*
\f}

where \f$\left[X_k\right]\f$ is the molar concentration of species \f$k\f$ and \f$\mathbf{V}_{\mathbf{k}}^*\f$ its molar diffusion velocity. The computation of the diffusion velocity \f$V_k\f$ is a central aspect of each transport model, and depends on the computed diffusion coefficients. 

#### Mixture-averaged transport model

The mixture-averaged model employs averaging rules on the pure species transport and binary diffusion coefficients to compute transport coefficients of species \f$\mathrm{k}\f$ into the mixture. The model offers a practical compromise between accuracy and computational expense.

For the computation of the mixture-averaged viscosity, the semi-empirical expression developed by Wilke and modified by Bird et al. is used

\f{align}{
\mu=\sum_{k=1}^N \frac{X_k \mu_k}{\sum_{k=1}^N X_j \Phi_{k j}},
\f}

where \f$\Phi_{k j}\f$ is defined by
\f{align}{
\Phi_{k j}=\frac{1}{\sqrt{8}}\left(1+\frac{W_k}{W_j}\right)^{-1 / 2}\left[1+\left(\frac{\mu_k}{\mu_j}\right)^{1 / 2}\left(\frac{W_j}{W_k}\right)^{1 / 4}\right]^2 .
\f}
For the mixture-averaged thermal conductivity, a similar expression is derived 

\f{align}{
\lambda=\frac{1}{2}\left(\sum_{k=1}^N X_k \lambda_k+\frac{1}{\sum_{k=1}^N X_k / \lambda_k}\right) .
\f}

For the formulation of the mixture-averaged diffusion coefficients, an equivalent diffusion coefficient for each species \f$k\f$ into the mixture of gases is computed. This mixture averaged diffusion coefficient is denoted by the subscript \f$m\f$. When computing the mixture-averaged diffusion coefficients, different expressions are available. An overview of the relations, as well as the exact expressions to compute the mixture averaged diffusion coefficients, is presented in the following table: 

| Flux | Diffusion Coefficient | Diffusion Velocity |
| :--- | :--- | :--- |
| \f$\mathbf{J}_k^*=-c D_{k m}^* \nabla X_k\f$ | \f$D_{k m}^*=\frac{1-X_k}{\sum_{j \neq k}^N X_j / \mathcal{D}_{k j}}\f$ | \f$\mathbf{V}_{\mathbf{k}}^*=-\frac{c}{X_k} D_{k m}^* \nabla X_k\f$ |
| \f$\mathbf{j}_k=-\rho D_{k m} \nabla Y_k\f$ | \f$\frac{1}{D_{k m}}=\sum_{j \neq k}^N \frac{X_j}{\mathcal{D}_{k j}}+\frac{X_k}{1-Y_k} \sum_{j \neq k}^N \frac{Y_j}{\mathcal{D}_{k j}}\f$ | \f$\mathbf{V}_k=-\frac{D_{k m}}{Y_k} \nabla Y_k\f$ |
| \f$\mathbf{j}_k=-\rho \frac{Y_k}{X_k} D_{k m}^{\prime} \nabla X_k\f$ | \f$D_{k m}^{\prime}=\frac{1-Y_k}{\sum_{j \neq k}^N X_j / \mathcal{D}_{k j}}\f$ | \f$\mathbf{V}_k=-\frac{1}{X_k} D_{k m}^{\prime} \nabla X_k\f$ |

The mixture-averaged expressions are first-order approximations, and unlike in the Stefan-Maxwell equations and in the multicomponent model, no constraint exists to ensure that the net species diffusion flux is equal to zero, e.g.

\f{align}{
\sum_{k=1}^N \rho Y_k V_k=0 .
\f}

#### Multicomponent transport model

The Chapman-Enskog theory provides the basis for the multicomponent model. This model is obtained after a rigorous approach to the derivation of the transport coefficients. The thermal diffusion coefficient, the thermal conductivity and the multicomponent diffusion coefficients are evaluated from the solution of a system of equations, which are defined by the matrix \f$L\f$. The matrix \f$L\f$ can be further decomposed into nine block submatrices. The equation system takes the form
\f{align}{
\left(\begin{array}{ccc}
L^{(00,00)} & L^{(00,10)} & 0 \\
L^{(10,00)} & L^{(10,10)} & L^{(10,01)} \\
0 & L^{(01,10)} & L^{(01,01)}
\end{array}\right)\left(\begin{array}{c}
a_{00}^1 \\
a_{10}^1 \\
a_{01}^1
\end{array}\right)=\left(\begin{array}{c}
0 \\
\mathbf{X} \\
\mathbf{X}
\end{array}\right) .
\f}

where the vector \f$\mathbf{X}\f$ is composed of the single mole fraction components \f$X_k\f$. The multicomponent diffusion coefficients are calculated from the inverse of the submatrix \f$L^{(00,00)}\f$ as

\f{align}{
D_{j k}=X_j \frac{16 T}{25 p m_k} \bar{W}\left(P_{j k}-P j j\right),
\f}

with

\f{align}{
(P)=\left(L^{(00,00)}\right)^{-1}
\f}

The thermal conductivity is found from the solution of the system of equations as
\f{align}{
\begin{aligned}
\lambda_{0, \text { trans }} & =-4 \sum_{k=1}^N X_k a_{k 10}^1 \\
\lambda_{0, \text { int }} & =-4 \sum_{k=1}^N X_k a_{k 01}^1 \\
\lambda_0 & =\lambda_{0, \text { trans }}+\lambda_{0, \text { int }}
\end{aligned}
\f}
and the thermal diffusion coefficient is computed as
\f{align}{
D_k^T=\frac{8 m_k X_k}{5 R} a_{k 00}^1 .
\f}
Similarly to the mixture-averaged model, the diffusion velocity can be computed through the equation
\f{align}{
V_k=\underbrace{\frac{1}{X_k \bar{W}} \sum_{j \neq k}^N W_j D_{k j} d_j}_{\text {Ficksian diffusion }}-\underbrace{\frac{D_k^T}{\rho Y_k} \frac{1}{T} \nabla T}_{\text {Soret diffusion }},
\f}
where \f$q_j\f$ is the diffusion driving force, which can be computed as
\f{align}{
d_j=\underbrace{\nabla X_j}_{\text {pecies gradient diffusion }}+\overbrace{\left(X_j-Y_j\right) \frac{1}{p} \nabla p}^{\text {Pressure gradient diffusion }} .
\f}


The necessary components of the \f$L\f$ matrix are

\f{align}{
\begin{aligned}
& L_{j k}^{00,00}=\frac{16 T}{25 p} \sum_{l=1}^K \frac{X_l}{m_j \mathcal{D}_{j l}}\left\{m_k X_k\left(1-\delta_{j l}\right)-m_j X_j\left(\delta_{j k}-\delta_{k l}\right)\right\}, \\
\hspace 2mm
& L_{j k}^{00,10}=\frac{8 T}{5 p} \sum_{l=1}^K X_k X_l\left(\delta_{j k}-\delta_{j l}\right) \frac{m_l\left(1.2 C_{k l}^*-1\right)}{\left(m_k+m_l\right) \mathcal{D}_{k l}}, \\
\hspace 2mm
& L_{j k}^{10,00}=L_{k j}^{00,10}, \\
\hspace 2mm
& L_{j k}^{01,00}=L_{k j}^{00,01}=0, \\
\hspace 2mm
& L_{j k}^{10,10}=\frac{16 T}{25 p} \sum_{l=1}^K \frac{m_j}{m_k} \frac{X_j X_l}{\left(m_j+m_l\right)^2 \mathcal{D}_{j l}}\left\{\left(\delta_{k l}-\delta_{j k}\right)\left[\frac{15}{2} m_k^2+\frac{25}{4} m_l^2-3 m_l^2 B_{j l}^*\right]\right. \\
& \left.-4 m_k m_l A_{j l}^*\left(\delta_{k l}+\delta_{j k}\right)\left[1+\frac{5}{3 \pi}\left(\frac{C_{j, \text { rot }}}{R \xi_{j l}}+\frac{C_{l, \text { rot }}}{R \xi_{l j}}\right)\right]\right\} \text {, } \\
\hspace 2mm
& L_{j j}^{10,10}=-\frac{16 m_j X_j^2}{R \mu_j}\left(1+\frac{10 C_{j, \mathrm{rot}}}{R \xi_{j j}}\right)-\frac{16 T}{25 p} \sum_{l \neq j}^K \frac{X_j X_l}{\left(m_j+m_l\right)^2 \mathcal{D}_{j l}} \\
\hspace 2mm
& \times\left\{\frac{15}{2} m_j^2+\frac{25}{4} m_l^2-3 m_l^2 B_{j l}^*+4 m_j m_l A_{j l}^*\left[1+\frac{5}{3 \pi}\left(\frac{C_{j, \mathrm{rot}}}{R \xi_{j l}}+\frac{C_{l, \mathrm{rot}}}{R \xi_{l j}}\right)\right]\right\}, \\
\hspace 2mm
& L_{j k}^{10,01}=\frac{32 T}{5 \pi p C_{k, \text { int }}} \sum_{l=1}^K \frac{m_k A_{k l}^*}{\left(m_k+m_l\right) \mathcal{D}_{k l}}\left(\delta_{j l}+\delta_{j k}\right) X_k X_l \frac{C_{k, \mathrm{rot}}}{R \xi_{k l}}, \\
\hspace 2mm
& L_{j j}^{10,01}=\frac{16}{3 \pi} \frac{m_j X_j^2}{\mu_j C_{j, \text { int }}} \frac{C_{j, \text { rot }}}{R \xi_{j j}}+\frac{32 T R}{5 \pi p C_{j, \text { int }}} \sum_{l \neq j}^K \frac{m_j A_{j l}^*}{\left(m_j+m_l\right) \mathcal{D}_{j l}} X_j X_l \frac{C_{j, \text { rot }}}{R \xi_{j l}}, \\ 
& L_{j k}^{01,10}=L_{k j}^{10,01}, \\ 
\hspace 2mm
& L_{j j}^{01,01}=-\frac{8 R^2}{\pi C_{j, \text { int }}^2} \frac{m_j X_j^2}{R \mu_j} \frac{C_{j, \text { rot }}}{R \xi_{j j}}-\frac{4 R T}{C_{j, \text { int }} p}\left\{\sum_{l=1}^K \frac{X_j X_l}{\left[\mathcal{D}_{j \mathrm{int}, l}\right]}+\sum_{l \neq j}^K \frac{12 X_j X_l}{5 \pi C_{j, \text { int }}} \frac{m_j}{m_l} \frac{A_{j l}^*}{\mathcal{D}_{j l}} \frac{C_{j, \mathrm{rot}}}{\xi_{j j}}\right\}, \\
\hspace 2mm
& L_{j k}^{01,01}=0(j \neq k) . \\ 
\end{aligned}
\f}


For a closer look at the computation and terms, please refer to [[Kee]]

[Kee]:https://onlinelibrary.wiley.com/doi/book/10.1002/9781119186304