# Gravitational Lensing Lab

MATLAB code for the gravitational **microlensing** experiment of the
undergraduate Physics Lab C course. It takes the light curve of a star that
was temporarily magnified by a compact object passing in front of it
(OGLE survey data), fits the point-source point-lens (Paczyński) model to
it, and estimates the physical parameters of the event, their
uncertainties, and the mass of the lens.

## Purpose

When a massive object (the *lens*) crosses the line of sight to a
background star (the *source*), its gravity bends the light and the star
appears brighter for days to months. For a point lens the flux, in units of
the un-lensed baseline flux, is

```
u(t)  = sqrt( u_min^2 + ((t - T0) / tau)^2 )
mu(u) = (u^2 + 2) / (u * sqrt(u^2 + 4))
F(t)  = 1 - f_bl + f_bl * mu(u(t))
```

| parameter | meaning |
|-----------|---------|
| `T0`      | time of closest approach (peak of the light curve) |
| `tau`     | Einstein-radius crossing time – sets the duration, and through it the lens mass |
| `uMin`    | impact parameter in Einstein radii – sets the peak magnification |
| `fbl`     | blending: fraction of the baseline light that really comes from the lensed star (1 = no blending) |

The code answers the questions of the lab, for any event:

1. **When and how bright is the peak?** – weighted least-squares parabola
   through the points around the maximum (`fitPeakParabola`).
2. **What are `tau` and `uMin`?** – χ² on a 2-D grid with `T0` fixed
   (`gridFit`).
3. **…and `T0`?** – 3-D grid.
4. **…and blending `fbl`?** – 4-D grid.
5. **How well are they determined?** – Δχ² contours for every pair of
   parameters (`plotChi2Contours`) and 1σ errors from Δχ² = 1.
6. **Cross-check** with a non-linear simplex fit (`fitPSPL`, `fminsearch`).
7. **What is the lens mass?** – from `tau`, assuming a transverse velocity
   and distances (`lensMass`).

## Usage

```matlab
main                       % analyse every data file listed in main.m
runtests('tests')          % unit tests (synthetic data with known answers)
```

Put the OGLE `phot.dat` files in `data/` (see `data/README.md`) and edit
the table at the top of `main.m`. With no data files present, `main` runs
the same analysis on a simulated event so you can see the output.

Single steps can be used on their own:

```matlab
addpath src
lc    = loadLightCurve("data/BLG002.dat");          % magnitudes -> normalised flux
guess = estimateInitialParams(lc);                  % rough T0, tau, uMin from the data
res   = gridFit(lc, guess, ["T0" "tau" "uMin" "fbl"]);
res.best, res.errLow, res.errHigh
plotChi2Contours(res);
plotLightCurve(lc, res.best);
lensMass(res.best.tau)                              % solar masses
```

Requires MATLAB R2021a or newer (no toolboxes).

## Layout

| file | role |
|------|------|
| `main.m` | driver: list of events → `analyzeEvent` for each |
| `src/analyzeEvent.m` | the whole pipeline for one light curve, with printout and figures |
| `src/loadLightCurve.m`, `lightCurveFromMagnitudes.m` | read a file, convert magnitude → baseline-normalised flux |
| `src/psplMagnification.m`, `uFromMagnification.m` | the model and its inverse |
| `src/estimateInitialParams.m` | starting values straight from the data |
| `src/polyfitWeighted.m`, `fitPeakParabola.m` | weighted polynomial fit with covariance; peak time and flux |
| `src/chi2Grid.m` | χ² on a grid of any 1–4 parameters |
| `src/gridFit.m` | self-refining grid search, best fit and 1σ errors |
| `src/gridMinimum.m`, `coordinateDescentMin.m` | N-dimensional minimum search (global / axis-by-axis) |
| `src/fitPSPL.m` | `fminsearch` fit with errors from the χ² curvature |
| `src/plotChi2Contours.m`, `plotLightCurve.m` | figures |
| `src/lensMass.m` | `tau` → lens mass |
| `src/makeSyntheticLightCurve.m` | simulated light curves for tests and the demo |
| `tests/testMicrolensing.m` | unit tests |
| `results/original/` | figures produced for the lab report with the original scripts |

## What changed relative to the original scripts

The first version of this code (still in the git history) had one script
per lab part and per number of parameters. It was merged into one general
toolbox, and these problems were fixed:

- **One N-dimensional implementation** instead of copies for 2, 3 and 4
  parameters (`chi2Grid`, `gridMinimum`, `coordinateDescentMin`,
  `plotChi2Contours` work for any subset of parameters).
- **No hard-coded events**: file paths (`~/Desktop/...`), the index window
  `493:533` of the peak, per-event search ranges and catalogue starting
  values are gone. Starting values are estimated from the data and the
  grid re-centres and rescales itself around the minimum.
- **Baseline normalisation** used the first data point (`flux(1)`), so one
  noisy measurement rescaled the whole light curve. It is now the median
  of all points far from the event.
- **Swapped indices**: `chi_mat(y, x, z)` was used to read the minimum of
  an array indexed `(x, y, z)`, so Δχ² was measured from the wrong
  element; in the 4-parameter fit the indices of `tau`/`uMin` and of
  `T0`/`fbl` were exchanged when converting back to parameter values.
- **Grid values vs. axes**: best-fit values were computed as
  `lower + index * delta` with `delta = range / N`, while the grid itself
  was `linspace(lower, upper, N)` – a different spacing. Values are now
  read from the grid vectors.
- **χ² accumulated across panels**: the six 2-D slices of the 4-parameter
  contour figure reused one matrix without resetting it.
- **`T0` search range** was `T0 * (1 ± 1e-6)` (an accident of the size of
  Julian dates) and the simplex start was `T0 * 0.99` (24 000 days away
  from the event). `T0` is now handled relative to the peak in units of
  `tau`.
- **Coordinate-descent minimum** can stop in a local minimum; the global
  minimum of the grid is used for the fit (the axis-by-axis routine is
  kept, generalised to N dimensions).
- **Uncertainties** are now computed (Δχ² = 1 on the profile χ², and from
  the Hessian for the simplex fit) rather than only read off the plots.
  Contours by default show the *profile* χ² (other parameters
  re-optimised), which gives the correct joint confidence regions; the
  old fixed-slice behaviour is available with `'Mode', 'slice'`.
- **Speed and memory**: the model is linear in `fbl`, so the 4-D grid
  costs the same as a 3-D one.
- Negative flux errors (sign of the derivative), an unfinished function
  (`findMaximum.m`), editor autosave files and duplicated copies were
  removed.

## References

- B. Paczyński, *Gravitational Microlensing in the Local Group*, ARA&A 34 (1996).
- R. Narayan & M. Bartelmann, *Lectures on Gravitational Lensing*.
- OGLE Early Warning System, <https://ogle.astrouw.edu.pl/ogle3/ews/ews.html>.
