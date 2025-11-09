# Levelset (LS) # {#mmLS}

The evolution of the level-set function 
\f$ \varphi(\mathbf{x},t) \f$ is described 
by the transport equation
\f{align}{
\frac{\partial \varphi}{\partial t} + v_{\Gamma} \cdot \nabla \varphi = 0 \quad, 
\f}

where \f$ \varphi \f$ represents the signed distance to a given interface.
The interface can be of a solid type, as for the wall- distance application 
or of an arbitrary shape with changing curvature, i.e., representing a 
flame front for combustion applications.

## Wall distance

For the wall distance application of a solid interface, 
the shape of the level-set is conserved and the level-set function is only 
transported in space.

## Flame front
For the modelling of the interaction between flame-front and flow field, special care has to be taken. In general, the formulation of the level-set equation permits an arbitrary, sufficiently smooth function \f$ \varphi(\mathbf{x},t) \f$ . However, solving the transport equation moves the zero level set \f$ \varphi_0 \f$ correctly, but may perturb the level-set function near \f$ \varphi_0 \f$, i.e., it may cause very large or small gradients. Various schemes have been developed to allow for the correct reinitialization of the problem and avoid the unphysical perturbation of the level-set function. To alleviate this difficulty, the arbitrary level-set function is replaced by a well behaved function and $\varphi$ is initialized into a signed distance function, which is the unique viscosity solution of the Eikonal equation
\f{align}{
|\nabla \varphi|=1
\f}

anchored at \f$\varphi_0\f$. However, once initialized into such a signed distance function, the level-set function \f$\varphi\f$ usually does not retain this property under the evolution of the level-set equation and needs to be reinitialized at regular time intervals. The most straightforward but inefficient reinitialization technique is to directly compute the minimum distance of each point from the zero level set, requiring a computational work of the order \f$\mathcal{O}\left(\mathcal{N}^2\right)\f$, with \f$\mathcal{N}\f$ being the number of cells. A more efficient and simpler approach is to use a partial differential equation to iteratively reinitialize the level-set function. Sussman et al. [[Sussman]] reformulate the Eikonal equation as an evolution equation in artificial time \f$\tau\f$

\f{align}{
\partial_\tau \varphi^\nu+S(\tilde{\varphi})\left(\left|\boldsymbol{\nabla} \varphi^\nu\right|-1\right)=0,
\f}

where the superscript \f$\nu\f$ denotes the discrete pseudo-time step. The quantity \f$S(\tilde{\varphi})\f$ is a smoothed sign function of the perturbed level-set function \f$\tilde{\varphi}=\tilde{\varphi}(\mathbf{x}, \tau=0)\f$ being defined as

\f{align}{
S(\tilde{\varphi})=\frac{\tilde{\varphi}}{\sqrt{\tilde{\varphi}^2+\epsilon^2}},
\f}

where \f$\epsilon\f$ is a smoothing parameter. Analytically, the modified equations yield for \f$\tau \rightarrow \infty\f$ the unique viscosity solution of the Eikonal equation correcting the perturbed level-set function \f$\tilde{\varphi}\f$ to become a signed distance function and keep the zero level set invariant because \f$S\left(\tilde{\varphi}_0\right)=0\f$.   

However, since in a discrete representation of \f$\tilde{\varphi}\f$ hardly any computational points coincide with \f$\tilde{\varphi}_0\f$, the location of the zero level set has to be defined by interpolating neighboring points, considerably displacing the zero level-set and thus may lead to substantial errors due to the reinitialization. 

In summary, a reinitialization method should have the following two properties:
* modify the level-set function such that \f$|\nabla \varphi|=1\f$ is satisfied,
* keep the \f$\varphi_0\f$ iso-surface invariant.

Both criteria can be fulfilled analytically, but one usually faces an overdetermined problem in the discrete version of the level-set function in multi-dimensional space such that a solution which exactly meets both criteria cannot be obtained.

Different formulations have been made to further circumvent this problem. The first formulation (**CR-1**) is based on the least-squares method which minimizes the unwanted displacement of the interface within the reinitialization. The overdetermined problem, which is solved in this first formulation of the reinitialization, is reduced to a determined problem in a second formulation such that the location of the interface is preserved within the reinitialization. The second formulation (**CR-2**) is derived by minimizing the number of constraints imposed on the reinitialization scheme in the first formulation. The formulations are modifications of the differential equation based methods introduced in [[Russo]] and modified in [[Sussman]] and shown above. Several other reinitialization methods are implemented as shown [here](@ref ugLS).

## References
G. Russo and P. Smereka. A remark on computing distance functions. J. Comput. Phys., 163:51–67, 2000. [https://doi.org/10.1006/jcph.2000.6553][Russo]

M. Sussman, P. Smereka, and S. Osher. A level set approach for computing solutions to incompressible two-phase flow. J. Comput. Phys., 114:146–159, 1994. [https://doi.org/10.1006/jcph.1994.1155][Sussman]

[Sussman]:https://doi.org/10.1006/jcph.1994.1155
[Russo]:https://doi.org/10.1006/jcph.2000.6553