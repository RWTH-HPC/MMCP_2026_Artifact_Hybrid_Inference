# Boundary conditions # {#mmBoundaryConditions}

At the domain boundaries the boundary conditions (BC) need to be applied to each face/wall defining the simulation domain, i.e. the walls surrounding the computational domain as well as solid bodies, like a sphere, airfoil, etc. inside the flow.
Boundary conditions can be subdivided into inlet boundary conditions, outlet boundary conditions, no-slip boundary conditions, constant pressure boundary conditions, axisymmetric boundary conditions, symmetric boundary conditions, and periodic or cyclic boundary conditions. [[Versteeg1995]]

Besides boundary conditions initial conditions need to be applied for transient simulations, defining initial values of the flow variables in the computationalg domain. 


## Dirichlet boundary condition
The Dirichlet boundary condition indicates the value of the solution function \f$\phi\f$ to the differential equation on the domain boundary \f$\partial \Gamma \f$ to match the known function \f$ \psi \f$.
\begin{eqnarray}  
\phi\vert_{\partial \Gamma} = \psi\vert_{\partial \Gamma}  
\end{eqnarray}  
Note that the function \f$ \psi \f$ might also depend on variables other than the surface coordinates, e.g. the time \f$t\f$.  
The no-slip condition for viscous fluids is an example for a Dirichlet boundary condition, where the velocity of the fluid normal to the wall vanishes and the tangential fluid velocity equals the wall velocity. For the slip condition only the velocity normal to the wall is set to zero whereas the velocity parallel to the wall is let free. 

## Neumann boundary condition 
Neumann boundary conditions force the normal derivative of the solution function \f$\phi\f$ to the differential equation on the domain boundary \f$\partial \Gamma \f$ to match the known function \f$ \psi \f$.
\begin{eqnarray} 
\fracpart{\phi}{n}\Bigg|_{\partial \Gamma} = \psi\vert_{\partial \Gamma}  
\end{eqnarray} 
Again \f$ \psi \f$ might be time-dependent. 
The Neumann boundary condition is used for application of symmetry planes or for the modeling of wall friction, if it is proportional to the strain rate. 


## Cauchy, Robin and Mixed boundary conditions
Those three boundary conditions are combinations of the Neumann and Dirichlet boundary condition. The Cauchy boundary condition specifies both the value and the normal derivative on the boundary, whereas the Robin boundary condition specifies a weighted average of Neumann and Dirichlet boundary conditions. In the mixed boundary condition the domain boundary is divided into different parts which need to fulfill either Neumann or Dirichlet boundary condition. 


## References
* Versteeg, H.K., and Malalasekera, W., An Introduction to Computational Fluid Dynamics - The Finite Volume Method., Longman Science & Technical., 1995, [Versteeg1995].

[Versteeg1995]: http://ftp.demec.ufpr.br/disciplinas/TM702/Versteeg_Malalasekera_2ed.pdf
