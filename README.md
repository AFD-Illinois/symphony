#*symphony*

##Overview

*symphony* calculates synchrotron emissivities and absorptivities, polarized in the Stokes basis {I, Q, U, V}, for any arbitrary gyrotropic momentum space distribution function.  Three distribution functions are built in: a relativistic thermal (Maxwell-Juettner) distribution, a nonthermal power-law distribution, and a kappa distribution.  

##Features

* `C` code to calculate synchrotron emissivities via the `j_nu()` function and absorptivities via the `alpha_nu()` function.
* `C` code to evaluate approximate fitting function values for the emissivity and absorptivity, via the `j_nu_fit()` and `alpha_nu_fit()` functions, respectively.
* CMake configure system, which helps during the build process to find all necessary libraries and files.
* `Python` interface for `j_nu()`, `alpha_nu()`, `j_nu_fit()`, and `alpha_nu_fit()`.
  * This combines the speed of `C` when evaluating emissivities and absorptivities with `Python`'s user-friendly syntax.  It also allows for interfacing with larger `Python` codes.

##How to use *symphony*

1. To download and build *symphony*:
  1. Clone *symphony* from github.  Navigate into the *symphony* folder, and note that it contains a folder named "src/"; create a folder named "build" and navigate into it.
  2. Type "cmake" followed by the location of the "src/" folder.  Altogether, this line should look something like: "cmake /location/to/symphony/src"
  3. Type "make"
  4. Try to run the *symphony* demo by typing "./symphony".  This should output three columns: "nu/nu_c", "j_nu()", and "j_nu_fit()".  Each column should have 8 rows, including the column titles mentioned above.  The three columns correspond to the frequency of emission divided by the cyclotron frequency, the calculated emissivity for some sample parameters, and the value of the fitting formula for those parameters.

2. To use *symphony*'s `Python` interface:
  1. Navigate to the "build/" folder created in step 1., above.  Open `Python` in the command line or by writing a ".py" file.
  2. Import `symphony` by typing "import symphonyPy".  This allows one to call 4 functions: `j_nu_py()`, `alpha_nu_py()`, `j_nu_fit_py()`, and `alpha_nu_fit_py()`.  The first two provide calculated values of the emissivity and absorptivity for the input parameters, and the latter two provide the corresponding approximate fitting formula results.
  3. To find the arguments of these functions, open up the docstring for the function you wish to use.  To do this (for example, for the function `j_nu_py()`), open up `Python` in the command line, type "import symphonyPy", hit enter, and type "symphonyPy.j_nu_py?".  This should display the docstring for `j_nu_py()`.  

* Arguments for `j_nu()` (in `C`) and `j_nu_py()` (`Python`): `j_nu(nu, magnetic_field, electron_density, observer_angle, distribution, polarization, theta_e, power_law_p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width)`.
 * Note: the arguments for `j_nu_py()`, `alpha_nu_py()`, `j_nu_fit_py()`, and `alpha_nu_fit_py()` are all identical.
 * The arguments are also all the same for the four `C` functions the above `Python` functions call.
* Sample call to `j_nu_py()`: `j_nu_py(230e9, 30, 1, 1.047, symphonyPy.MAXWELL_JUETTNER, symphonyPy.STOKES_I, 10, 2.5, 1, 1000, 1e10, 3.5, 10)`
 * Note: All parameters with units are in CGS.
