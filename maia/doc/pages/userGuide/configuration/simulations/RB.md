# Rigid bodies (RB) # {#ugRB}

[TOC]

The rigid bodies (RB) solver is used to track finite sized bodies for the simulation of fully-resolved particle laden problems.
Typically, this solver is coupled with a flow solver of your choice to simulate fluid particle interaction.
Here, the particles are tracked in a Lagrangian frame of reference. Despite the Lagrangian nature of this solver, it is based on an Cartesian grid for the following reasons:

- Spatial partitioning of the particles
- Efficient coupling with other Eulerian solvers
- Parallelization

The solver offers a number of different shapes as rigid bodies. For the solver itself, shape and size only matter for the parallelization, as the solution itself only describes the center-of-mass values of each body. To couple the solver to any flow solver, the actual surface of each body can be provided by this solver.

## Running a simulation

To enable this solver, you must define a solver with the type
```toml
solverType = "MAIA_RIGID_BODIES"
```

For a single settling particle an example configuration looks like this:

```toml
initialBodyCenters  = [2.8, 0.0, 0.0]  # Position x,y,z in STL units
bodyRadius          = 0.5              # Radius in STL units
bodyDensityRatio    = 1.167            # Adapt if Re is changed
uRef                = 1.053            # Adapt if Re is changed
integrateRotation   = false           
forcedMotion        = false            # Integrate motion equation OR use prescribed motion
gravity             = [9.81, 0.0, 0.0] # Gravity vector
noEmbeddedBodies    = 1                # Number of bodies used
```

### Body initialization

If multiple bodies are used, their initial position can be given by

```toml
initialBodyCenters = [x_b1, y_b1, z_b1, ..., x_bn, y_bn, z_bn]
```
Alternatively, a box seeding can be used to spawn multiple bodies in an equidistant manner. The following configuration snippet spawns a total of 3x3x3=27 bodies.

```toml
bodyCenterInitMethod  = "BOX_SEED"
seedBoxMin            = [-3.0, -3.0, -3.0]
seedBoxMax            = [3.0, 3.0, 3.0]
bodiesPerEdge         = 3
maxNoBodies           = 30
```
