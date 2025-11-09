# LPT - X # {#nmLptX}

The Lagrangian particle tracking (LPT) solver requires flow field information from a flow solver.
To support non-congruent grids, the LPT solver itself stores the flow field information corresponding to the used grid.
If two-way models are applied, the source terms for mass, energy, and momentum are calculated in the LPT and need to be transferred to the flow solver.
The transfer of flow field information and source terms as well as their potential interpolation between the different grids is done by the coupler.

\dot
digraph D {
  rankdir=LR;
  node [shape=record];
  fs [label = "Flow solver"];
  c [label = "Coupler"];
  subgraph cluster_lpt {
    label = "LPT";

    node [shape=record];
    rank = same;

    cells [label = "Cells"];
    parts [label = "Particles"];

    cells -> parts [dir = "both"];
  }

  fs -> c [label = "flow field"];
  c -> cells [label = "flow field"];
  cells -> c [label = "source terms"];
  c -> fs [label = "source terms"];
}
\enddot

Two main different approaches for the timesteping of the LPT solver and the flow solver, 
in this example of a finite volume type, exist. 
A simple sucessive time stepping is possible, which means that one solver time step is 
executed after the other. If a blocking exchange occoures in the solver sub-steps, 
such as a particle or source term exchange in the LPT solver, this time stepping can be exessively inefficient. 
As an alternative, an interleafed execution is possible, in which the different solver communication 
can be hidden behind the execution of the sub-step of the other solver.
The coupling of a 3-step flow solver and a LPT time-stepping split into an equal number of sub-steps is 
displayed as an example for the two different apporaches in the figure below. 

![](LPTX-interleafed.svg) { width=60% }
