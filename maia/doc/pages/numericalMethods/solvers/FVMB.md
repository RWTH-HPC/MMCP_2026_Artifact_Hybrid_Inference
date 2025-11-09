# Finite volume - moving boundary (FVMB) # {#nmFVMB}

In the following, the numerical method developed for the simulation of flows with freely moving
boundaries is described based on the solution scheme published in [[Schneiders1, Schneiders2]. 
The method enables a sharp resolution of the embedded boundaries 
and strictly conserves mass, momentum, and energy. 
A new explicit Runge-Kutta scheme (PC-RK) is used for these applications as a predictor-corrector
type reformulation of a popular class of Runge-Kutta methods which substantially reduces the
computational effort for moving boundaries tracking and subsequent solver reinitialization.
The solver stability and accuracy is not impaired.


## Runge-Kutta Scheme 
The conservative quantities \f$\mv{Q} \f$ are integrated from time level \f$ t^n \f$ to 
\f$t^{n+1} = t^{n + \Delta t}\f$ by this new explicit Runge-Kutta scheme. 
This new formulation is a modification of the widely-used class
of low-storage schemes originally proposed by Van Der Houwen [87, 88] and later established
in computational aerodynamics by Jameson [93].
These multi-stage Runge-Kutta schemes (MS-RK), however, have a significant
drawback for problems in which the residual operator itself changes over time such as in moving
boundary problems. Then, to sustain a time-consistent formulation and the accuracy of the
scheme, the residual operator \f$ R(t; \cdot)\f$ has to be reconstructed at each intermediate time level
defined by a Runge-Kutta substage, i.e., at \f$t = t^n + \alpha_{k−1} \Delta t, k = 1, . . . , s\f$. 
Since this usually involves surface tracking and remeshing, 
the latter corresponds to the computation of the cut-cell geometries, as well as a subsequent reinitialization of the solver, e.g., re-computing the
least-squares system, it renders the scheme inefficient when this overhead is on the
order of the time required for evaluating \f$R(t; \mv{Q})\f$ and applying the MS-RK.
As a remedy, we propose a new class of explicit schemes in which the intermediate time levels
are eliminated while strictly retaining the stability properties and accuracy the MS_RK. Using
the same RK-coefficients as for the MS_RK and letting \f$\alpha_s ≡ 1\f$ 
for time consistency, the single-time-level
counterpart to the MS-RK is
\f{align}{
(\mv{Q}V)^{(0)} &= (\mv{Q}V)^n , \\
(\mv{Q}V)^{(1)} &= (\mv{Q}V)^{(0)} − \Delta t \, R(t^n ; \mv{Q}^{(0)}), \\
(\mv{Q}V)^{(k)} &= (\mv{Q}V)^{(0)} − \Delta t [(1 − \alpha_{k−1}) R(t^n ; \mv{Q}^{(0)}) + \alpha_{k−1} R(t^{n+1}; \mv{Q}^{(k−1)})] \quad , k = 2, . . . , s \\
(\mv{Q}V)^{n+1} &= (\mv{Q}V )^{(s)}.
\f}
The first stage, corresponds to a forward-Euler prediction step while the succeeding
stages increase the stability of the scheme by applying convex combinations of the residuals
\f$R(t^n; ·)\f$ and \f$R(t^{n+1}; ·)\f$ using the step size \f$\Delta t\f$. 
In this predictor-corrector Runge-Kutta scheme (PC-RK) all intermediate values, 
\f$ \mv{Q}^{(0)} , \mv{Q}^{(0)} , . . . , \f$ are approximations to the solution at time level
\f$t^{n+1}\f$. More precisely, each substep of the PC-RK contains the same residual operators \f$R(t^n ; ·)\f$
and \f$R(t^{n+1} ; ·)\f$, i.e., it is not varied between the substages. Assuming that the initial residual
\f$R(t^n; \mv{Q}^n )\f$ has been computed and stored at the previous time step, the PC-RK scheme involves
only residual evaluations at time level \f$t^{n+1}\f$ . Hence, only a single construction of the residual
operator, i.e., a single remeshing and reinitialization has to be performed per Runge-Kutta cycle.
This reduces the overall costs of constructing \f$R\f$ by a factor of \f$s\f$ compared to the MS-RK scheme. The costs of storing \f$R(t^n ; \mv{Q}^n)\f$ and the additional operations of the PC-RK are
usually negligible in comparison.
Let \f$t_{\text{init}}\f$ denote the time required to construct the residual and \f$t_{\text{exec}}\f$ the time to execute the
new scheme, i.e., s-times evaluation of the residual and application of the RC-RK. The speedup
of the PC-RK over the MS-RK can be estimated by \f$1 + (s − 1)\sigma \f$, where \f$\sigma = t ini t /(t_{\text{init}} + t_{\text{exec}} ) \f$ is
the computational overhead for the construction of the residual in the new scheme. The magnitude of \f$\sigma \f$ depends mainly on the ratio of cut cells to regular cells. For the cases investigated
in this paper, \f$\sigma\f$ varies between 0.1 and 0.38. For \f$s = 5\f$, the new scheme therefore is about
1.4 − 2.5 times faster, showing the effort for the moving-boundary tracking to be drastically reduced.
Moreover, the coupling of the new scheme to other solvers, e.g., for structural motion,
Lagrangian particle tracking, a level-set method, or heat conduction, is remarkably simplified
since boundary conditions are only required at the time level \f$t^{n+1}\f$, 
i.e., the interfacial states do not have to be resolved at the intermediate time levels by the coupled solver. 
