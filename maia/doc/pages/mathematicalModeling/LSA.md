# Linear scalar advection # {#mmLSA}

The linear scalar advection equation 

\f{equation}{
  \frac{\partial \phi} {\partial t} + \nabla \cdot (\mv{u} \phi) = 0
\f}

describes the advection of a scalar \f$\phi\f$ with a constant advection velocity vector \f$\mv{u}\f$. The initial condition for scalar \f$\phi \f$ is given by \f$ \phi(t=0) = \phi_0 \f$.  
Due to its simplicity this equation is computationally very cheap and can be used for testing and developing features. Further it serves as an evaluation tool for the Discontinuous Galerkin (DG) operator. 

  In the m-AIA solver, the linear scalar advection equation is solved in the class MAIADdSysEqnLinearScalarAdv with MAIADgSysEqnLinearScalarAdv::m_advectionVelocity being an n-dimensional input parameter prescribed in the "properties_run.toml" script used for running the simulation. 
