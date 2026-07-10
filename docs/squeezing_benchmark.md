# CW parametric-squeezing benchmark (flat pump)

An exact, closed-form check of the undepleted squeezing dynamics, used by the
`"CW parametric squeezing (flat pump)"` testset in `test/runtests.jl`.

## Model

`Linearized{Squeezing{QD3WM}}` is undepleted-pump vacuum squeezing: a fixed
classical second-harmonic (SH) pump drives the fundamental (FH) signal. In
normalized units, with signal GVD $\beta_2 = s$, coupling $\varepsilon$, and a
**flat** (translationally invariant) pump $\beta_x = \beta_0$,

$$
\mathrm{i}\,\partial_z \hat\Phi_x = -\tfrac{s}{2}\,\partial_x^2 \hat\Phi_x
      + \varepsilon\beta_0\,\hat\Phi_x^\dagger .
$$

Because the pump is flat, the dynamics is **diagonal in frequency**: a flat
pump carries zero momentum, so it couples only the signal pair
$(\omega, -\omega)$. Each pair is an independent two-mode squeezer with
generator

$$
M(\omega) = \begin{pmatrix} D & \gamma \\ -\gamma & -D \end{pmatrix},
\qquad D(\omega) = \tfrac{s}{2}\,\omega^2, \qquad \gamma = \varepsilon\beta_0 ,
$$

and $M^2 = (D^2 - \gamma^2)\,I$.

## Exact result

Starting from vacuum and propagating a length $Z$, mode $\omega$ is transformed
by the two-mode-squeezer Bogoliubov coefficients

$$
A(\omega) = \cosh(\nu Z) - \mathrm{i}\,\tfrac{D}{\nu}\sinh(\nu Z),
\qquad
B(\omega) = -\mathrm{i}\,\tfrac{\gamma}{\nu}\sinh(\nu Z),
\qquad
\nu = \sqrt{\gamma^2 - D^2},
$$

where a complex $\nu$ gives amplification ($\sinh/\cosh$) for $D < \gamma$ and
oscillatory sidelobes ($\sin/\cos$) for $D > \gamma$. The two physical
observables are the photon-number (gain) spectrum and the entangled-pair
correlation:

$$
n(\omega) = \langle \hat a_\omega^\dagger \hat a_\omega \rangle = |B|^2
\quad\text{(number sector)},
\qquad
m(\omega) = \langle \hat a_\omega \hat a_{-\omega} \rangle = \mathrm{i}\,A\,B
\quad\text{(pair sector)},
$$

with $n(0) = \sinh^2(\gamma Z)$. The spectrum $n(\omega)$ is amplified near
$\omega = 0$ and develops **$\sin^2$ sidelobes with exact zeros** once
$D(\omega) > \gamma$ — strongly non-monotonic. The $\mathrm{i}$ in $m(\omega)$
is the package's anomalous-correlator phase convention (a single global factor,
verified empirically).

On the discrete grid each Fourier mode is *exactly* such a two-mode squeezer,
so this holds mode-by-mode (not merely in a continuum limit). Hence:

- the **normal** covariance $\Pi = \langle \delta\hat\Phi^\dagger \delta\hat\Phi \rangle$
  (`a⁺a`) is **diagonal** in the Fourier basis (each mode marginally thermal —
  the entangled partners have zero number-sector coherence,
  $\langle \hat a_\omega^\dagger \hat a_{-\omega} \rangle = 0$), with diagonal
  $n(\omega)$;
- the **anomalous** covariance $\Sigma = \langle \delta\hat\Phi\,\delta\hat\Phi \rangle$
  (`aa`) is **anti-diagonal** ($\langle \hat a_\omega \hat a_{\omega'} \rangle
  = W\Sigma W^{\mathsf T}$ is nonzero only at $\omega' = -\omega$) — this is where
  the $\omega \leftrightarrow -\omega$ entanglement lives — with anti-diagonal
  $m(\omega)$.

Diagonality of $\Pi$ is thus the *fingerprint* of two-mode squeezing, not
evidence against entanglement; the entanglement is carried entirely by $\Sigma$.

## What the test checks (all exact to the RK4IP floor)

Parameters $N=128,\ L=20,\ s=1,\ \varepsilon=1,\ \beta_0=1,\ Z=0.5,\ N_\text{steps}=500$:

Number sector ($\Pi$):

1. $\Pi$ diagonal in the Fourier basis (off-diagonal fraction $< 10^{-10}$) —
   translational invariance.
2. $n(0) = \sinh^2(\gamma Z)$ — the DC amplification anchor (fixes
   $\gamma = \varepsilon\beta_0$).
3. $n(\omega)$ matches $|B|^2$, *including the $\sin^2$ sidelobes*, to
   $< 10^{-5}$ (observed $\sim 8\times10^{-8}$ at 500 steps, converging
   $\sim N_\text{steps}^{-4}$).

Pair sector ($\Sigma$):

4. $\Sigma$ anti-diagonal in the Fourier basis (off-anti-diagonal fraction
   $< 10^{-10}$) — each $\omega$ entangles with $-\omega$ and nothing else.
5. the full complex $m(\omega)$ matches $\mathrm{i}\,A\,B$, *phase included*, to
   $< 10^{-5}$.

Any error in the coupling constant, dispersion factor / $2\pi$ convention,
parametric sign, vacuum term, or the $\Sigma/\Pi$ dispersion coupling moves
$n(\omega)$ or $m(\omega)$ by $O(1)$ and is caught.
