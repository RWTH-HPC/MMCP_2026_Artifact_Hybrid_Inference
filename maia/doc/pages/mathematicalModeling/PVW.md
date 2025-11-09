# Structural mechanics # {#mmFC}

[TOC]

In the following, the principle of virtual work and the material laws of linear elastic solids
are briefly described. For further information, the reader is referred to [].

## Principle of virtual work (PVW) # {#mmFCPVW}
The FCM derives from the principle of virtual work, which is based on the equilibrium of internal
stresses \f$ \sigma \f$ and external volume forces \f$ f \f$ in a solid body:
\f{equation}{
  div(\mv{\mv{\sigma}}) + \mv{f} = 0 \f}
The virtual work \f$ \delta W \f$ is defined as the product of real forces with virtual
displacements and vis versa. Multiplying the equilibrium with a virtual displacement
\f$ \delta \mv{u} \f$ and integrating over the body volume yields
\f{equation}{
  \delta W = \int\limits_V (\mv{f} + div(\mv{\mv{\sigma}})) \cdot \delta \mv{u} \, dV \f}
Using the product rule the equation can be reformulated yielding
\f{equation}{
  \delta W = \int\limits_V \mv{f} \cdot \delta \mv{u} + div(\mv{\mv{\sigma}}\delta \mv{u}) - \mv{\mv{\sigma}} : \mv{\nabla} \delta \mv{u} 
  \, dV \f}
where the operator \f$ : \f$ denotes the Frobenius scalar product of the stress tensor
\f$ \mv{\mv{\sigma}} \f$ and the tensor \f$ \mv{\nabla} \delta \mv{u} = grad(\delta\mv{u}) \f$.
Using the divergence theorem, the equation can be rewritten as
\f{equation}{
  \delta W = \int\limits_V \mv{f} \cdot \delta \mv{u} \, dV + \int\limits_A \mv{t} \delta \mv{u} \, dA - 
    \int\limits_V \mv{\mv{\sigma}} : \mv{\nabla} \delta \mv{u} \, dV \f}
where \f$ t \f$ is the traction vector on the body surface. For bodies in equilibrium, the
virtual work is \f$ \delta W = 0 \f$. Hence, the weak form of the principle of virtual work
is given by
\f{equation}{
  \int\limits_V \mv{f} \cdot \delta \mv{u} \, dV + \int\limits_A \mv{t} \delta \mv{u} \, dA = \int\limits_V \mv{\mv{\sigma}} :
    \mv{\nabla} \delta \mv{u} \, dV \f}

## Material laws # {#mmFCMaterial}
So far, only linear elastic materials are considered by the FC solver. Furthermore, it is limited
to small strains. For such materials, the relation between displacements, strains, and stresses
is given by
\f{equation}{
  \varepsilon_x = \frac{\partial u_x}{\partial x} \f}
and
\f{equation}{
  \sigma_x = E \varepsilon_x \f}
in 1D. In 2D and 3D the strains and stresses are described using the matrices
\f{equation}{
  \mv{\mv{\varepsilon}} = \left(\begin{array} . 
  \varepsilon_{x} & \varepsilon_{xy} & \varepsilon_{xz} \\
  \varepsilon_{yx} & \varepsilon_y & \varepsilon_{yz} \\ 
  \varepsilon_{zx} & \varepsilon_{zy} & \varepsilon_z \end{array} \right) \f}
and
\f{equation}{
  \mv{\mv{\sigma}} = \left(\begin{array} . 
  \sigma_{x} & \sigma_{xy} & \sigma_{xz} \\
  \sigma_{yx} & \sigma_y & \sigma_{yz} \\ 
  \sigma_{zx} & \sigma_{zy} & \sigma_z \end{array} \right) \quad .\f} 
The notation used in the code uses the symmetry of the matrices
\f$ \varepsilon_{ij} = \varepsilon_{ji} \f$, writing them as vectors.
\f{equation}{
  \mv{\varepsilon} = \left(\begin{array} .
      \varepsilon_x \\
      \varepsilon_y \\
      \varepsilon_z \\
      \varepsilon_{xy} \\
      \varepsilon_{xz} \\
      \varepsilon_{yz} \end{array} \right) \f}
\f{equation}{
  \mv{\sigma} = \left(\begin{array} .
      \sigma_x \\
      \sigma_y \\
      \sigma_z \\
      \sigma_{xy} \\
      \sigma_{xz} \\
      \sigma_{yz} \end{array} \right) \f}
Hence, the relation between strains and stresses is given by
\f{equation}{
  \mv{\sigma} = \mv{\mv{C}} \mv{\varepsilon} \f}
with
\f{equation}
  \mv{\mv{C}} = \left(\begin{array} .
    c_1 & c_2 & c_2 & 0 & 0 & 0 \\
    c_2 & c_1 & c_2 & 0 & 0 & 0 \\
    c_2 & c_2 & c_1 & 0 & 0 & 0 \\
    0 & 0 & 0 & c_3 & 0 & 0 \\
    0 & 0 & 0 & 0 & c_3 & 0 \\
    0 & 0 & 0 & 0 & 0 & c_3 \end{array} \right) \f}
The coefficients of matrix \f$ C \f$ derive from Hook's law and are given by
\f{equation}{
  c_1 = \frac{E (1 - \nu)}{(1 - \nu) (1 - 2\nu)} \f}
\f{equation}{
  c_2 = \frac{\nu}{(1 - \nu) (1 - 2\nu)} \f}
and
\f{equation}{
  c_3 = \frac{E}{1 + \nu} \f}
