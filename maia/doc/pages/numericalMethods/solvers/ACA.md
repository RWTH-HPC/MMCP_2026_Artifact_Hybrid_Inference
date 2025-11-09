# Acoustic Analogy (ACA) # {#nmACA}

## Non-dimensionalization # {#nmACA_nonDim}
In the ACA solver the variables are made non-dimensional with a infinity state,
see [Tab.1](@ref nmACA_nonDim_refVarsTable). The non-dimensionalization of the
relevant variables is given in [Tab.2](@ref nmACA_nonDim_nonDimVarTable).

<table>
<caption id="nmACA_nonDim_refVarsTable">Tab.1: Reference variables</caption>
<tr><td>Free-stream speed of sound  </td> <td>\f$ c_\infty \f$                          </td></tr>
<tr><td>Free-stream density         </td> <td>\f$ \rho_\infty \f$                       </td></tr>
<tr><td>Reference length            </td> <td>\f$ L_{ref} = 1.0 \f$ (in STL file units) </td></tr>
</table>

<table>
<caption id="nmACA_nonDim_nonDimVarTable">Tab.2: Non-dimensionalization in ACA solver</caption>
<tr><th>Variable                    </th> <th>Equation                                                           </th></tr>
<tr><td>Spatial coordinates         </td> <td>\f$  \mv{x}^* = \mv{x}  \cdot \frac{1}{L}                       \f$</td></tr>
<tr><td>Perturbed velocity          </td> <td>\f$ \mv{u'}^* = \mv{u'} \cdot \frac{1}{c_\infty}                \f$</td></tr>
<tr><td>Perturbed density           </td> <td>\f$   \rho'^* =  \rho'  \cdot \frac{1}{\rho_\infty}             \f$</td></tr>
<tr><td>Perturbed pressure          </td> <td>\f$      p'^* =     p'  \cdot \frac{1}{\rho_\infty c_\infty^2}  \f$</td></tr>
<tr><td>Time                        </td> <td>\f$       t^* =     t   \cdot \frac{c_\infty}{L_{ref}}          \f$</td></tr>
</table>

## Fourier transformation (FT) # {#nmACA_fourierTransformation}
Solving the FWH equation numerically is computational more efficient in
frequency domain [[Lockard2000]]. The Fourier transformation is a mathematical
tool to perform this transformation and is defined in the following.

A signal in time domain \f$ q(t) \f$ is transformed into frequency domain
\f$ Q(\omega) \f$ using the Fourier transformation \f$ \mathcal{F}(\cdot) \f$,
\f{equation*}{
  \mathcal{F}(q(t)) := \int_{-\infty}^{\infty} q(t) e^{-\iu\omega t} dt = Q(\omega)
  \,,
\f}
where \f$ t \f$, \f$ f \f$, and \f$ \omega = 2\pi f \f$ are the time, frequency,
and angular frequency, respectively.<br>
The inverse Fourier transformation \f$ \mathcal{F^{-1}}(\cdot) \f$, i.e., the
transformation from frequency into time domain, is defined as
\f{equation*}{
  \mathcal{F}^{-1}(Q(\omega)) := \frac{1}{2\pi}\int_{-\infty}^{\infty} Q(\omega) e^{\iu\omega t} d\omega = q(t)
  \,.
\f}
## Discrete Fourier transformation (DFT)
The [Fourier transformation](@ref nmACA_fourierTransformation) is defined for
functional signal over an infinite sampling duration. However, practical
application provide discrete and finite samplings making the FT impractical.

The discretized Fourier transformation (DFT) is defined for a finite duration
\f$ T \f$ with \f$ N \f$ equidistant samples spaced with
\f$ \Delta t = \frac{T}{N} \f$. Defining the discrete frequency \f$ f_n \f$ and
time \f$ t_k \f$
\f{align*}{
  f_n &:= \frac{n}{N\dt} , \, \text{for} \quad
  n \in \mathbb{N} \land n\in \left[ -\frac{N}{2}, \frac{N}{2} \right] \,,\\
  %
  t_k &:= k \Delta t , \, \text{for} \quad
  k \in \mathbb{N} \land n\in [0, N-1] \,,\\
\f}
and introducing the abbrevations
\f{align*}{
  Q_n &:= Q(f_n) \\
  q_k &:= q(t_k)
\f}
the DFT reads
\f{equation*}{
  Q_n
    = \sum_{k_0}^{N-1} q_k e^{\iu 2\pi f_n t_k} \dt 
    = \dt \sum_{k_0}^{N-1} q_k e^{\iu 2\pi \frac{k~n}{N}} \,.
\f}
The inverse DFT is
\f{equation*}{
  q_k
    = \frac{1}{N} \sum_{n=0}^{N-1} Q_n e^{-\iu 2\pi \frac{k~n}{N}}
\f}

### Fast Fourier transformation (FFT)
The implementation of a DFT is straightforward but its computation is cumbersome
for bigger problems. Hence, for number of samples corresponding to a power of
two ( \f$ N=2^i \,, \text{for} \, i\ \in \mathbb{N} \f$ ) the Fast Fourier
Transformation (FFT) algorithm [Press1992] is a more efficient alternative.

## Windowing
While in the FT an infinite signal length is required the DFT allows for a
finite but periodic signal as input. If the provided sample is not consisting of
multiple oscillation periods distortions are introduced - often called spectral
leakage.
This disturbance can be reduced by multiplying the signal with window function,
which zeros the signal towards its start and end. While a rectangle window
describes the raw signal the following functions are common
\f{align*}{
  \text{0. Rectangle       : } w_k &= 1 \\
  \text{1. Hanning         : } w_k &= \frac{1}{2} \left( 1 - \cos \left( \frac{2\pi k}{N} \right) \right) \\
  \text{2. Hamming         : } w_k &= 0.54 - 0.45 \cos \left( \frac{2\pi k}{N} \right) \\
  \text{3. Modified Hanning: } w_k &= 
  \begin{cases}
    1 \,, &\text{for} \, k \in \left[ \frac{N}{8}, \frac{7N}{8} \right] \\
    \frac{1}{2} \left( 1 - \cos \left( \frac{8\pi k}{N} \right) \right) \,, &\text{else}
  \end{cases}
\f}

Applying the window function, i.e., multiplying it with the signal, reduces the
energy contained in the signal. Energy conservation is ensured by scaling the
data after the Fourier transformation by the factor \f$ 1/E_T \f$
\f{equation*}{
  E_T = \sqrt{ \frac{1}{N} \sum_{k=0}^{N-1} |w_k|^2} \,.
\f}

## Implementation

\dot
digraph G {
  rankdir=TB;

	subgraph cluster_0 {
		label = "Generate input (not ACA solver)"
		color="#6BAD00";
		node [shape=rectangle];
    CFDCAA[label="CFD/CAA"];
    Analy[label="Analytically genereated data"];
    GenerateStl[label="Generate surface STL"];
    SampledData[label="Sampled data with respect to STL"];
    ObserverPFromFile[label="Generate observer points from file"];

    GenerateStl -> CFDCAA;
    GenerateStl -> Analy;
    CFDCAA -> SampledData;
    Analy -> SampledData;
	}

	subgraph cluster_1 {
		label = "ACA Solver";
		color="#006DB2";
		node [shape=rectangle];
    ObserverPFromAnly[label="Generate observer points from analytical"];
    ChangeNonDim[label="Change non-dimenzionalization"];
    CalcSource[label="Calc source terms (Q,F)"];
    Window[label="Windowing"];
    FFT[label="FFT/DFT"];
    Integration[label="Calc integrals"];

    ChangeNonDim -> CalcSource;
    ObserverPFromAnly -> CalcSource;
    subgraph cluster_3 {
      label = "";
      CalcSource -> Window;
      Window -> FFT;
      FFT -> Integration;
    }

    subgraph cluster_2 {
      label = "Post processing";
      color="#BD0006";
      node [shape=rectangle];
      pabs[label="p'_abs"];
      prms[label="p'_rms"];
      spl[label="SPL"];
      oaspl[label="OASPL"];
    }
    Integration -> pabs;
    Integration -> prms;
    Integration -> spl;
    Integration -> oaspl;
	}


  ObserverPFromFile -> CalcSource;
  SampledData -> ChangeNonDim;

}
\enddot

## References
* Lockard, D. P., “An efficient, two-dimensional implementation of the Ffowcs Williams and Hawkings equation,” J. Sound Vibr.,
Vol. 229, No. 4, 2000, pp. 897–911. [doi.org/10.1006/jsvi.1999.2522][Lockard2000].
* W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, "Numerical Recipes
in C", Cambridge University Press, 1992.
[Lockard2000]: https://doi.org/10.1006/jsvi.1999.2522
