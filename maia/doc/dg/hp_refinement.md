Hp-refinement in the DG block
=============================

How does it work?
-----------------

### Uniform refinement

*   Occurs if all cells/elements have the same level and the same polynomial
    degree
*   Preparations for solver run:
    *   For each leaf cell (i.e. cell without children), a DG element is created
    *   Between two elements that share a face, a DG surface is created
*   During solver run:
    *   Values from the interior of the elements are interpolated to the element
        faces and saved to the Gauss nodes at the surfaces
    *   At each surface, a numeric flux is calculated using values from both the
        left and right elements
    *   The unique numeric flux is used to calculate the contributions of the
        surface integral to the time derivative in both neighboring elements
*   Special cases:
    *   MPI boundaries: no change in the algorithm, only need to exchange
        interpolated values on surface befure flux calculation
    *   Domain boundaries: surface has only one neighboring element; either
        outer state or flux are set by boundary condition


### P-refinement

*   Occurs if the polynomial degrees of the left and right element
    of a surface are different
*   No changes in preparations
*   During solver run:
    *   After interpolation, values at surface are at different locations from
        left and right element
    *   Before flux calculation, values from lower polynomial degree are
        interpolated to Gauss nodes of higher polynomial degree (forward
        projection)
    *   After flux calculation, fluxes are interpolated back on side with lower
        polynomial degree (reverse projection) before integral calculation


### H-refinement

*   Occurs if there is a Cartesian level difference between adjacent cells
*   Preparations for solver run:
    *   One surface is created for each *refined* element (i.e. at most `2` in
        2D, at most `4` in 3D)
*   During solver run:
    *   After interpolation, values at refined surfaces are at different
        locations
    *   Before flux calculation, values from coarse element are interpolated to
        Gauss nodes of refined element (forward projection)
    *   After flux calculation, fluxes are interpolated back from the refined
        surfaces to the the Gauss nodes of the coarse element (reverse
        projection) and summed up


### Hp-refinement

*   Occurs if h- and p-refinement appear at the same surface
*   No changes in preparations to normal h- or p-refinement
*   During solver run:
    *   Forward projection first for h-, then for p-refinement
    *   Reverse projection first for h-, then for p-refinement
    *   General principle: flux calculations at the highest level of
        hp-refinement to avoid loss in accuracy (i.e. all values are
        interpolated to the refined surfaces at the highest participating
        polynomial degree)


### Adaptive p-refinement

*   Occurs at predefined intervals before a time step calculation (see property
    `adaptiveInterval`)
*   Polynomial degree of each element can be increased or reduced by one, or
    kept constant
*   On change in polynomial degree, node coordinates have to be re-calculated
    and solutions interpolated to the new Gauss node locations



How to use it?
--------------

*   P-refinement: currently only hardcoded in DG block (`initPrefinement()`)
*   H-refinement: use grid generator to control local cell size through boundary
    refinement, refinement patches etc.
*   Adaptive p-refinement: set property `adaptiveRef` to one of supported
    methods (currently `DG_ADAPTIVE_TEST` and `DG_ADAPTIVE_GRADIENT`)


How is it implemented?
----------------------

### P-refinement
*   No changes to data structures
*   All elements/surfaces are allocated with space for `maxPolyDeg`
*   Algorithms must use local polynomial degree (instead of `initPolyDeg` or
    just the polynomial degree of the first element)
*   Result of forward projection is stored in temporary location and then
    copied back to surface
*   Result of reverse projection is stored in temporary location and then
    used for surface integral calculation

### H-refinement
*   Changes to data structures:
    *   Surfaces have a new field `fineCellId` to store the *cell* id of the
        refined element in case of h-refinement
    *   Additional collector `helements` stores additional surfaces for
        coarse elements in case of h-refinement
*   Before forward projection, results of interpolation are copied to all
    additional surfaces for coarse element
*   Separate loops for coarse and fine/unrefined elements (for coarse element,
    loops go over `helements`)
*   `fineCellId` used to determine surface location with respect to coarse
    element

### Adaptive p-refinement:
*   tbd.
