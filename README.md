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
1. Clone *symphony* from github.  Navigate into the *symphony* folder, and note that it contains a folder named "src/"; create a folder named "build" and navigate into it.
2. Type "cmake" followed by the location of the "src/" folder.  This should be something like: "cmake /location/to/symphony/src"
3. Type "make"
4. Try to run the *symphony* demo by typing "./symphony".  This should output three columns: "nu/nu_c", "j_nu()", and "j_nu_fit()".  Each column should have 8 rows, including the column titles mentioned above.  The three columns correspond to the frequency of emission divided by the cyclotron frequency, the calculated emissivity for some sample parameters, and the value of the fitting formula for those parameters.
