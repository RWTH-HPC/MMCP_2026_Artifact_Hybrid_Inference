# Finite Cell (FC) # {#nmFC}

[TOC]

This section outlines the underlining methods and algorithms of the FC solver.
This includes the discretiztion, the p-refinement, the subcell integration, and
the formulation of essential and natural boundary conditions.

## Discretiztion
The weak form of principle of virtual work ([PVW](@ref mmFCPVW)) is discretized by
using a shape function to map the distribution of the displacement to the whole cell.
The shape function is based on the Lagrangian polynominals. The Lagrangian polynominal
interpolates the the displacement between discrete nodes located inside the cell
also referred to as element. The number of nodes per element is defined by the
p-refinement, which is explained in the next section. The minimum number of nodes
is 2 in each space direction. Hence, in 3D, 8 nodes (one at each element vertex)
exist per element. The Lagrangian polynominal reads
\f{equation}
  N_{ijk}(x, y, z) = \prod_{0 \leq l \leq (p+2); l \neq i} \frac{x - x_i}{x_l - x_i}
  \prod_{0 \leq l \leq (p+2); l \neq j} \frac{y - y_j}{y_l - y_j}
  \prod_{0 \leq l \leq (p+2); l \neq k} \frac{z - z_k}{z_l - z_k}
  \f}
The subscribts \f$ i\f$, \f$ j\f$, and \f$ k \f$ denote the node position in 
\f$ x \f$, \f$ y \f$, and \f$ z \f$ direction. That is, \f$ N_{0,0,0} \f$ denotes the 
zeroth node position in \f$ x \f$, \f$ y \f$, and \f$ z \f$ direction. \f$ N_{1,2,1} \f$
denotes the first node position in \f$ x \f$, the second node position in \f$ y \f$
and the first node position in \f$ z \f$ direction. The separte shape functions are
assembled in the matrix \f$ \mv{\mv{N}} \f$. To obtain the strain-interpolation matrix
\f$ \mv{\mv{B}} \f$, the derivative of matrix \f$ \mv{\mv{N}} \f$ is calculated.
The linear system of equations solved by the FC solver reads
\f{equation}
  \mv{\mv{K}} \mv{u} = \mv{F} \f}
where \f$ \mv{\mv{K}} \f$ is the gobal stiffness matrix, \f$ \mv{u} \f$ is the global
displacement vector, and \f$ \mv{F} \f$ is the global force vector. The global vectors
and matrices are assembled from the element matrices/vectors. The element stiffness
matrix is defined by
\f{equation}
  \mv{\mv{K}} = \int\limits_{V_e} \mv{\mv{B}}^T \mv{\mv{C}} \mv{\mv{B}} det|\mv{\mv{X}}| \, dV_e \f}
where \f$ V_e \f$ is the element volume and \f$ det|\mv{\mv{X}}| \f$ is the determinant of the
Jacobian matrix.

## P-refinement # {#nmFCpRef}
When using p-refinement, the number of nodes per element is increased. That is,
increasing the p-refinement by \f$ 1 \f$ increases the number of nodes per Cartesian
direction by \f$ 1 \f$. The total number of nodes per element is given by
\f$ (p + 2)(p + 2)(p + 2) = p^3 + 6p^2 + 12p + 8 \f$. The distribution of the 
points follows the distribution of the Lobatto-Points or an equidistant distribution.

@warning Although the p-refinement is implemented until p = 10, it is tested so far only
for upto p = 3.

## Subcell integration # {#nmFCsubCell}
To increase the accuracy at boundaries without increasing the degree of freedom,
a subcell integration can be applied at curved boundaries. That is, the boundary
cell is subdivided into \f$ 4 \f$ (2D) or \f$ 8 \f$ (3D) cells. For each subcell
it is checked, if the subcell has a cut with the geometry. If yes, the cell is
further subdivided. If not, the integration is performed on the cell using 
\f$ \alpha \ll 1 \f$ if the cell is located outside and \f$ \alpha = 1 \f$ if
the cell is inside the geometry. The subdivision is performed until the maximum
subcell depth is reached. On the maximum depth, the inside-outside check is conducted
for each integration point and the \f$ \alpha \f$ is changed depending on the 
check.

@note Although the degree of freedome is not increased due to the subcell integration,
the computational cost are increased. For example, a subcell depth of 10 results in 
upto \f$ 10^9 \f$ subcells per boundary cell.

## Boundary conditions # {#nmFCBC}
In the following, the two types of boundaries are shown.

### Essential boundary conditions
At essential boundaries the displacement is set. This is done by solving the
following integral and adding the resulting matrix to the element stiffness
matrix.
\f{equation}
  \mv{\mv{K}} = \mv{\mv{K}} + \int\limits_A \mv{\mv{N}}^T k \mv{\mv{N}} det|\mv{\mv{X}}| \, dA \f}
Here, \f$ k \f$ is a factor with \f$ k \gg 1 \f$ and \f$ A \f$ is the boundary
surface, which is approximated by the triangles of the STL file. Essential boundary conditions
have the following numbers, which can be set in the `geometry.toml` file:
```
BC.0 = 8010

BC.0 = 8011

BC.0 = 8012
```

@note So far, only a displacement of \f$u = 0\f$ can be prescribed.

### Natural boundary conditions
At natural boundaries the surface traction is set. This is done by solving the
following integral and adding the resulting vector to the element force vector.
\f{equation}
  \mv{F} = \mv{F} + \int\limits_A \mv{\mv{N}}^T \mv{t} \, det|\mv{\mv{X}}| \, dA \f}
Here, \f$ t \f$ denotes the surface traction and \f$ A \f$ is the boundary
surface, which is approximated by the triangles of the STL file. Natural boundary conditions
have the following numbers, which can be set in the `geometry.toml` file:
```
BC.0 = 8030

BC.0 = 8031

BC.0 = 8032

BC.0 = 3035
```

