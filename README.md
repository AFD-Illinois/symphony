# symphony

## Overview
This code, *symphony*, calculates synchrotron emissivition and absorption coefficients polarized in the Stokes basis {I, Q, U, V}.  It does this explicitly for three momentum-space electron distributions: a relativistic thermal (Maxwell-Jüttner) distribution, a nonthermal power-law distribution, and a kappa distribution.  *symphony* is designed to work with any arbitrary distribution function, though, and the process to add new distribution functions should be straightforward and is outlined below.

## Features
* Code to calculate synchrotron emission and absorption coefficients polarized in each of the four Stokes parameters {I, Q, U, V}, for the relativistic thermal, power-law, and kappa distribution functions
* Fitting functions for each polarization/distribution function
* Straightforward functions j_nu(nu, B, n_e, theta) and alpha_nu(nu, B, n_e, theta) can be used to calculate the emissivity or absorptivity for a given set of parameters
  * Also included are functions j_nu_fit(nu, B, n_e, theta) and alpha_nu_fit(nu, B, n_e, theta), which produce values from the corresponding approximate fitting formulae, mentioned above, and are useful in cases where computation time is a constraint
