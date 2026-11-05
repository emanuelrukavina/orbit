# libOrbit

Analyze measured accelerations of two celestial bodies and checks whether an additional, unseen mass could explain the difference between their observed and predicted motion.

You provide the bodies' masses, positions, measured accelerations, known external accelerations, and measurement uncertainties. The library runs the analysis internally and returns one `Result` with a possible object's estimated mass and position, fit statistics, and an evidence classification.

**Current scope:** one observation time, one positive point mass, and a candidate position on the line segment between the two observed bodies. libOrbit does not reconstruct an orbit or prove that an object exists.

> Originally built this for my KSP mod, but since KSP 2 sucks, imma just publish the library here for safekeeping in case I ever need it again.

> UPDATE: Fuck..... https://www.ign.com/articles/take-two-shutters-kerbal-space-program-2-studio-amid-layoffs

## How it works

The library predicts how bodies A and B should accelerate under their mutual gravity and any known external influences. It subtracts these predictions from the measured accelerations. The remaining differences are called *residual accelerations*.

It then searches for a mass C between A and B whose gravitational pull best explains both residuals. Finally, it compares that fit against the model with no hidden mass, using simulated noise to estimate how surprising the improvement is.

## The math

Newton's law of gravitation is

$$
F = G\frac{m_1m_2}{r^2},
$$

where $G = 6.67430\times10^{-11}\ \mathrm{m^3\,kg^{-1}\,s^{-2}}$.

For positions $\mathbf r_A$ and $\mathbf r_B$, the acceleration of A caused by B is

$$
\mathbf a_{A\leftarrow B}=Gm_B\frac{\mathbf r_B-\mathbf r_A}{\lVert\mathbf r_B-\mathbf r_A\rVert^3}.
$$

The acceleration of B caused by A is calculated in the same way, with A and B exchanged. The library adds the supplied known external accelerations to obtain the predicted accelerations $\mathbf a_A^{(0)}$ and $\mathbf a_B^{(0)}$.

The unexplained parts of the observations are

$$
\Delta\mathbf a_A=\mathbf a_A^{\mathrm{obs}}-\mathbf a_A^{(0)},
\qquad
\Delta\mathbf a_B=\mathbf a_B^{\mathrm{obs}}-\mathbf a_B^{(0)}.
$$

Let $D=\lVert\mathbf r_B-\mathbf r_A\rVert$ and let $\hat{\mathbf n}=(\mathbf r_B-\mathbf r_A)/D$ point from A toward B. A candidate body C is placed at

$$
\mathbf r_C=\mathbf r_A+\lambda(\mathbf r_B-\mathbf r_A),
\qquad 0<\lambda<1.
$$

For example, $\lambda=0.5$ is the midpoint. A positive mass $m_C$ at that position would cause

$$
\mathbf g_A=\frac{Gm_C}{\lambda^2D^2}\hat{\mathbf n},
\qquad
\mathbf g_B=-\frac{Gm_C}{(1-\lambda)^2D^2}\hat{\mathbf n}.
$$

In the ideal case of exact, perfectly aligned residual accelerations, define

$$
\alpha=\Delta\mathbf a_A\cdot\hat{\mathbf n},
\qquad
\beta=-\Delta\mathbf a_B\cdot\hat{\mathbf n}.
$$

When $\alpha>0$, $\beta>0$, and both residuals have no perpendicular component, the inferred location and mass are

$$
\lambda=\frac{\sqrt{\beta}}{\sqrt{\alpha}+\sqrt{\beta}},
\qquad
m_C=\frac{D^2\alpha\beta}{G(\sqrt{\alpha}+\sqrt{\beta})^2}.
$$

Real observations contain noise, so libOrbit fits the mass and location instead of requiring an exact match. With acceleration uncertainties $\sigma_A$ and $\sigma_B$, it minimizes

$$
\chi^2(m_C,\lambda)=
\frac{\lVert\Delta\mathbf a_A-\mathbf g_A\rVert^2}{\sigma_A^2}
+
\frac{\lVert\Delta\mathbf a_B-\mathbf g_B\rVert^2}{\sigma_B^2}.
$$

For a fixed location, define the acceleration produced per unit mass:

$$
\mathbf s_A=\frac{G}{\lambda^2D^2}\hat{\mathbf n},
\qquad
\mathbf s_B=-\frac{G}{(1-\lambda)^2D^2}\hat{\mathbf n}.
$$

Then the best nonnegative mass at that location is

$$
\hat m_C(\lambda)=\max\left(0,\frac{P(\lambda)}{Q(\lambda)}\right),
$$

where

$$
P=\frac{\mathbf s_A\cdot\Delta\mathbf a_A}{\sigma_A^2}
+\frac{\mathbf s_B\cdot\Delta\mathbf a_B}{\sigma_B^2},
\qquad
Q=\frac{\mathbf s_A\cdot\mathbf s_A}{\sigma_A^2}
+\frac{\mathbf s_B\cdot\mathbf s_B}{\sigma_B^2}.
$$

The library searches the configured interval of $\lambda$ values and keeps the best fit.

To assess the fit, it compares the no-hidden-body error $\chi_0^2$ with the best hidden-body error $\chi_1^2$:

$$
\Delta\chi^2=\chi_0^2-\chi_1^2.
$$

A larger improvement is not automatically a discovery. libOrbit generates synthetic Gaussian-noise observations under the no-hidden-body model, repeats the full location search, and estimates a bootstrap p-value:

$$
p=\frac{1+N_{\mathrm{simulated}\ \Delta\chi^2\geq\mathrm{observed}\ \Delta\chi^2}}{B+1},
$$

where $B$ is the number of simulations. With bootstrap enabled, a positive-mass fit and $p\leq0.05$ are classified as `candidate_supported`. This is evidence **under the supplied model**, not proof of a real object. With `bootstrap_samples = 0`, the classification is `inconclusive`.

At a fixed location, the nominal mass uncertainty is $\sigma_m=1/\sqrt{Q}$, and the nominal 5σ mass sensitivity is $5/\sqrt{Q}$. These are **conditional on location**, not globally calibrated detection thresholds.

> The exact reconstruction equations explain the underlying physics. The implementation uses a constrained weighted least-squares fit, scanning candidate positions on a fixed grid and refining the best grid candidate numerically.

## Build

Requires a C++17 compiler and CMake.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/orbit_example
```

The project builds a static library. The public header is `include/liborbit.hpp`; implementation files and private headers are in `src/`.

## Use

Include the single public header and call `analyze()`:

```cpp
#include <liborbit.hpp>

int main() {
    liborbit::Input input;

    // Set input.a and input.b using the types in liborbit.hpp.
    // Provide masses, positions, measured accelerations,
    // known external accelerations, and uncertainties.

    liborbit::Result result = liborbit::analyze(input);

    // Inspect result for the estimated candidate and classification.
}
```

See `examples/main.cpp` for a complete example with actual values. All units are SI: kilograms, meters, and meters per second squared. Positions and accelerations must use the same inertial reference frame. The supplied external accelerations must **exclude** the mutual gravity of A and B, which the library calculates itself.

## Limitations

libOrbit currently works with a single observation epoch and only searches for one positive point mass between the observed bodies. It does not estimate velocity, fit a time-series orbit, or search general three-dimensional space. Known masses and positions are treated as exact; acceleration errors are assumed independent and Gaussian. Missing forces or incorrect uncertainty estimates can produce misleading results.

`candidate_supported` means the model improves the fit at the chosen bootstrap threshold; `no_significant_evidence` does not rule out all unseen objects. If the fitted mass is zero, its reported position has no physical meaning. This is a learning and experimentation library, not a validated astronomical detection instrument.