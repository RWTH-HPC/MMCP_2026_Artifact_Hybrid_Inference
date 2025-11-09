# Coupling # {#nmCoupling}

All solvers operate on a hierarchical Cartesian grid in which cells are organized in
an octree structure with parent, child, and neighbor relationships [[Hartmann2008a]].
The joint domain decomposition based on this underlying Cartesian grid
ensures an efficient spatial coupling between the different numerical 
methods. Solver cells contained in the same volume of the 
three-dimensional domain are assigned to the same process
allowing in-memory exchange of the coupling terms.
Individual spatial solver constraints for the
different physical systems are taken into account during the
mesh generation and solution adaptive grid refinement, 
where cells are tagged according to their use by the solver.
When different domain sizes are used, 
single solvers can be deactivated for certain subdomains, 
thus not participating in the solver communication.
<!---
%![](figures/coupling_SolverAffil.svg)
-->
An example of the spatial decomposition and solver affiliation for a
Direct Injection application is presented in the figure above.
The number of possible combinations of solver affiliations for a single cell 
is reduced by physical and application driven solver restrains, i.e.,
LPT and FV cells are both reduced to the fluid domain.
The LPT solver is additionally deactivated in the intake 
and exhaust ports due to the closed valves shortly
after start of injection.
Note that level-set cells are required 
outside the fluid domain for the resolution of the zero level-set contour.

@subpage nmDgX  
@subpage nmLptX  
@subpage nmEE  
@subpage nmZonal  
