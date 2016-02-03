# symphony

## Overview
This code, *symphony*, calculates synchrotron emission and absorption coefficients polarized in the Stokes basis {I, Q, U, V}.  It does this explicitly for three momentum-space electron distributions: a relativistic thermal (Maxwell-Jüttner) distribution, a nonthermal power-law distribution, and a kappa distribution.  *symphony* is designed to work with any arbitrary distribution function, though, and the process to add new distribution functions should be straightforward and is outlined below.

## Features
* Code to calculate synchrotron emission and absorption coefficients polarized in each of the four Stokes parameters {I, Q, U, V}, for the relativistic thermal, power-law, and kappa distribution functions
* Fitting functions for each polarization/distribution function
* Straightforward functions `j_nu(nu, B, n_e, theta)` and `alpha_nu(nu, B, n_e, theta)` can be used to calculate the emissivity or absorptivity for a given set of parameters
  * Also included are functions `j_nu_fit(nu, B, n_e, theta)` and `alpha_nu_fit(nu, B, n_e, theta)`, which produce values from the corresponding approximate fitting formulae and can be useful in cases where computation time is a constraint.

## How to use *symphony*

1. To directly use *symphony* to calculate emission/absorption coefficients, or to test the provided fitting formulae:
  1. Select the distribution function in file **symphony.h**.  Do this by setting the `#define DISTRIBUTION_FUNCTION` to `THERMAL` (for the relativistic thermal/Maxwell-Jüttner), `POWER_LAW` for the power-law, or `KAPPA_DIST` for the kappa distribution.
  2. Select the Stokes parameter in file **symphony.h**.  Do this by setting the `#define POL_MODE` to `STOKES_I`, `STOKES_Q`, `STOKES_U`, or `STOKES_V`
  3. Set parameters `B` (magnetic field strength), `n_e` (electron number density), `obs_angle` (observer angle, measured with respect to the magnetic field) in **main.c**
  4. `make` the files
  5. Running *symphony* at this point will output three columns: `nu/nu_c`, `ans`, and `fit`.  The first column outputs the frequency divided by the cyclotron frequency (nu_c = eB/2pimc); the second is the emissivity or absorptivity at that frequency for the set parameters; the third is the value of the corresponding fitting function for those parameters.

2. To call `j_nu()`, `alpha_nu()`, `j_nu_fit()`, or `alpha_nu_fit()` from another C code:
  1. Copy *symphony*'s files into the same directory as your code
  2. `#include 'symphony.h'` in your code
  3. [NEED TO FIX THIS] Select the distribution function and Stokes parameter as in steps 1.1 and 1.2 above in **symphony.h**
  4. Use the functions as you wish. Each has the same arguments.  For example, `j_nu(nu, B, n_e, theta)`
 
3. To include a new electron distribution function:
  1. Add the new distribution function as a function only of the Lorentz factor, gamma, to the file **calc.c**.  Thus the function should look like: `double new_dist(double gamma)`
  2. In the file **symphony.h** under the comment `\\choose distribution function`, add `#define NEW_DIST (3)` and change `#define DISTRIBUTION_FUNCTION (x)` to `#define DISTRIBUTION_FUNCTION (NEW_DIST)`
  3. Calculate (by hand) the differential in equation 13 of Pandya, Zhang, Chandra, and Gammie
  4. Immediately above `#elif DISTRIBUTION_FUNCTION == POWER_LAW` in the function `differential_of_f(double gamma, double nu)` in **calc.c**, add `#if DISTRIBUTION_FUNCTION == NEW_DIST` and enter the differential calculated in step 2 immediately below
  5. `make` *symphony*
  6. Try to run *symphony* with the new electron distribution function.  Debug as necessary.
