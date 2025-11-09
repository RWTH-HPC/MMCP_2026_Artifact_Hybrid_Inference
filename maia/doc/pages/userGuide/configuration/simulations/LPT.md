# Lagrangian particle tracking (LPT) # {#ugLPT}

[TOC]

The Lagrangian particle tracking (LPT) solver is used to track a finite number of point particles in a Lagrangian frame of reference
within a Eulerian phase.
Applications of the LPT solver are versatile and are ranging from spray injection including
particle breakup and evaporation over to particle-laden flows of small non-spherical biomass particles.

In the following an introduction on how to set up a simulation using the LPT solver is given.

## Running a simulation

Within the simulation framework m-AIA, the LPT solver is always coupled to one of the flow solvers.
You are therefore invited to also take a look at the user guide of the corresponding flow solver
(i.e. [FV](@ref ugFV) or [LB](@ref ugLB)).
The coupling of LPT and flow solver is described in more detail in the numerical methods
(c.f. [Coupling](@ref nmCoupling)).

In order to use the LPT solver, you must define a solver with the type

```toml
solverType = "MAIA_PARTICLE"
```

Simulations using the LPT solver be conducted using dimensionalized or non-dimensionalized quantities (**recommended**).
The latter is activated by setting the property

```toml
nonDimensionaliseLPT = true
```

### Particle initialization

In the general case, all particles are initialized with a predefined density and temperature, in the
non-dimensionalized case the specified values must be non-dimensionalized with the reference density and
temperature.
The particle density and temperature can be set as follows:

```toml
particleDensity     = 400   # in the non-dimensional case this is the ratio between particle and ref. density
particleTemperature = 1     # in the non-dimensional case this is the ratio between particle and ref. temperature
```

The initialization of particle size, position and velocity depends on the method chosen by the
`particleInitializationMethod`.
You can automatically initialize a given number of particles per cell.
One or more diameters can be specified, by default it is assumed the particle particles are equally distributed among
the given diameters. In addition, the initial particle velocity can be specified.
This initialization method can be used by defining the following properties:

```toml
particleInitializationMethod = 0
particleInitNoPartperCell    = xx                                 # set here the number of particles per cell
particleDiameter             = [diameter1, diameter2]
particleInitVelocity         = [velocityX, velocityY, velocityZ]
```

Another possibility is to initialize particles according to an input text file, which is especially useful for
testing purposes. For this purpose, set the initialization method to `particleInitializationMethod=4` and create a
file called `part.txt` which contains in each row the following initial particle information:

```
diameter density x-coordinate y-coordinate z-coordinate x-velocity y-velocity z-velocity
```

Further information on the implemented particle initialization methods can be found
[here](@ref particleInitializationMethod).

### Particle motion

In the LPT solver the particle movement is determined by integration of a given motion equation using a
[Predictor-Corrector scheme](@ref nmLPT).
Depending on the selected motion equation, additional analytical or empirical correction terms may be used.
An overview of the existing motion equations and drag laws are given [here](@ref motionEquation) and
[here](@ref particleDrag).

As a starting point for your simulation the non-dimensional motion equation for inertial spherical particles using the
Schiller-Naumann drag correction is recommended.
This can be done by setting

```toml
motionEquation = 1   # Non-dimensional motion equation
particleDrag   = 2   # Schiller-Naumann drag term
```

Gravity effects are accounted for by specifying the Froude number

```toml
Frm = [0.0, 0.0, 0.0]
```

@note The Reynolds number specified by the property ```Re``` and used by the LPT solver corresponds to the Reynolds number
based on the stagnation point (\f$Re_0\f$) rather than the infinity state (\f$Re_\infty\f$).
<br> The coupled flow solver might have another convention, i.e. for the finite volume solver this property refers to
the infinity state.

### Coupling and evaporation

The LPT solver also supports **two-way coupling**, i.e., the modelling of fluid-particle interaction.
In @maia, coupling mechanisms for momentum, mass and energy heat are implemented.
In addition, evaporation can be considered.
The added source terms can be redistributed to neighboring cells by activating the source term redistribution and
specifying the number of redistributions layers.
You can specify the corresponding properties as follows:

```toml
particleMomentumCoupling = true
particleMassCoupling     = false   # specifies if mass coupling should be considered
particleHeatCoupling     = true    # specifies if heat coupling should be considered
particleCouplingRedist   = true    # specifies if the source term should be redistributed to neighboring cells
particleNoRedistLayer    = 2       # specifies the number of layers for the redistribution
```

### Particle-Wall and Particle-Particle collisions

The LPT solver supports particle-Wall as well as Particle-Particle collisions (**four-way coupling**).
Particle-Wall collisions can be activated by specifying

```toml
particleWallCollisions = true
```

Particle-Particle can be specified using the `particleCollisions`.
Further information on the implemented modes is given [here](@ref particleCollisions).

### Ellipsoidal particles

The LPT solver also offers Lagrangian point-particle models for particles with an ellipsoidal shape.
The ellipsoidal particles can coexist next to spherical particles and have their own motion equations and drag laws
as for the ellipsoidal point-particle model additional degrees of freedom must be considered due to the particle
orientation.
Ellipsoidal particles are activated using the property `particleIncludeEllipsoids`.
A rotational symmetric ellipsoid is well-defined by specifying its semi minor axis and a corresponding aspect ratio.

```toml
particleIncludeEllipsoids      = true
particleEllipsoidsemiMinorAxis = [semiMinorAxis1, semiMinorAxis2, ...]
particleEllipsoidaspectRatio   = [aspectRatio1, aspectRatio2, ...]         # aspectRatio1 refers to semiMinorAxis1
```

## Adaptive mesh refinement

If you want to learn how to enable AMR for you simulation, read [here](@ref ugAMR). The list of all available sensors for the LPT solver can be found @subpage sensorsLPT "here".

### Example

The code snippet below shows a minimal example of the properties required by the LPT solver to conduct a simulation
coupled with the Finite Volume solver.
The LPT solver solves the non-dimensional motion equation for spherical particles extended with the Schiller-Naumann
drag correction for the Lagrangian phase.
In each cell of the domain, a particle with a diameter of \f$0.01\f$ is spawned with \f$0\f$ velocity.
A gravitational force, defined through the `Frm` property, acts on the particles in the z-direction.
The definition of the Froude number is given in the [mathematical modeling](@ref mmLPT).

<details>
<summary>Example of the properties related to the LPT solver</summary>

```toml
noSolvers                    = 2
multiSolverGrid              = true
solverOrder_0                = [0,1]
solverOrder_1                = [1,0]
recipeMaxNoSteps             = 2
couplerOrder_0               = [1,1]
adaptationOrder              = [0,1]
noCouplers                    = 1
solversToCouple_0            = [0,1]

solvertype.0                 = "MAIA_FINITE_VOLUME"     # flow solver i.e. here finite volume
solvertype.1                 = "MAIA_PARTICLE"          # LPT Solver
couplerType_0                = "COUPLER_FV_PARTICLE"    # coupling of flow solver (i.e. FV) and LPT

motionEquation               = 1                        # non-dimensional motion equation
particleDrag                 = 2                        # Schiller-Naumann drag
nonDimensionaliseLPT         = true
particleInitializationMethod = 0
particleInitNoPartperCell    = 1.0
particleDiameter             = [0.1]
particleInitVelocity         = [0.0, 0.0, 0.0]
Frm                          = [0.0, 0.0, 0.32]
```

</details>

### Misc.

Depending on your setup, the following additional properties could be of interest:

<table>
  <tr>
    <th>Property</th>
    <th>Default</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>particleSizeLimit</td>
    <td>1E-12</td>
    <td>Particles smaller than the defined size are removed.</td>
  </tr>
  <tr>
    <td>maxNoParticles</td>
    <td>1E6</td>
    <td>Defines the maximum number of particles per rank.</td>
  </tr>
  <tr>
    <td>particleExhangeBufferSize</td>
    <td>1000</td>
    <td>Defines the size of the buffer used to exchange particles.</td>
  </tr>
  <tr>
    <td>particleRespawn</td>
    <td>false</td>
    <td>
      Flag to respawn particles that left the domain. <br>
      Particles are respawned on a plane defined by the property particleRespawnPlane
    </td>
  </tr>
</table>
