# *symphony*

## Overview

*symphony* calculates synchrotron emissivities and absorptivities, polarized in the Stokes basis {I, Q, U, V}, for any arbitrary gyrotropic momentum space distribution function. As of 2018, *symphony* also computes Faraday rotation coefficients in the Stokes basis and the full relativistic magnetized plasma susceptibility tensor, provided the distribution function is isotropic. Three distribution functions are built in: a relativistic thermal (Maxwell-Juettner) distribution, a nonthermal power-law distribution, and a kappa distribution.

*symphony* was first described in [Pandya et al., 2016 ApJ 822 34](http://dx.doi.org/10.3847/0004-637X/822/1/34), and extended to compute Faraday rotation coefficients in [Pandya et al., 2018](https://arxiv.org/abs/1810.05530). If you use this code in an academic context, cite these papers.

## Features

* `C` code to calculate synchrotron emissivities via the `j_nu()` function and absorptivities via the `alpha_nu()` function.
* `C` code to calculate Faraday rotation coefficients via the `rho_nu()` function
* `C` code to evaluate approximate fitting function values for the emissivity and absorptivity, via the `j_nu_fit()` and `alpha_nu_fit()` functions, respectively.
* CMake configure system, which helps during the build process to find all necessary libraries and files.
* `Python` interface for `j_nu()`, `alpha_nu()`, `j_nu_fit()`, and `alpha_nu_fit()`.
* `Python` interface to `rho_nu()` and spline fits
  * This combines the speed of `C` when evaluating emissivities and absorptivities with `Python`'s user-friendly syntax.  It also allows for interfacing with larger `Python` codes.

## How to use *symphony*

### To download and build *symphony*:
 1. Clone *symphony* from github.
 2. Unzip the "kernel_samples_datafiles.zip" file in place within the root "symphony" folder.
 3. Edit the "susceptibilityPy.pyx" file within `src/susceptibility_tensor/` so that the `main_directory` variable (on line 288) points to the absolute path of the "kernel_samples_datafiles/" folder that was produced in step 2.
 4. Create a new folder named "build" within the root directory and navigate into it.
 5. Type `cmake` followed by the location of the "src/" folder in the root directory. Altogether, this line should look something like: `cmake /location/to/symphony/src`. You can add the argument `-DCMAKE_INSTALL_PREFIX=/name/of/dir` to set the name of the directory to install to. You can use the relative paths here if you like.
 6. Type `make`.
 7. Optionally, run `make install` to install the library and Python module onto your system.
 8. In the "build" folder, navigate into the newly created "susceptibility_tensor" folder. Make an empty file called "\__init\__.py" and save it.

### To use *symphony*'s `Python` interface:
 1. Navigate to the "build/" folder created in step 1., above.  Open `Python` in the command line or by writing a ".py" file.
 2. Import *symphony* by typing `import symphonyPy`.  
  * This allows one to call the functions: `j_nu_py()`, `alpha_nu_py()`, `j_nu_fit_py()`, `alpha_nu_fit_py()`, `rho_nu_py()`, and spline fit functions `alpha_{I, Q, V}_spline()`, `rho_{Q, V}_spline()`, and `chi_ij()`.  
  * The first two provide calculated values of the emissivity and absorptivity for the input parameters, and the second two provide the corresponding approximate fitting formula results.
  * The function `rho_nu_py()` computes Faraday rotation coefficients using the *symphony* integrator; this is slow, and should only be used as a check of the spline fits.
 3. The arguments of these functions can be found by accessing the associated docstrings.  This can be done in the `Python` command line using the following: 
```
import symphonyPy
symphonyPy.j_nu_py?
```

### Arguments for Emissivity, Absorptivity, and Faraday Rotation Coefficient Functions
* Arguments for all emissivity, absorptivity, and Faraday rotation coefficient functions in both `C` and `Python` take nearly the same arguments and output a double.  The only difference is that the C version has an `error_message` parameter for handling evaluation errors, and the absorptivity coefficients include a switch to choose the integration method. The arguments are: 
```
j_nu(nu, magnetic_field, electron_density, observer_angle, distribution, polarization, 
     theta_e, power_law_p, gamma_min, gamma_max, gamma_cutoff, kappa, kappa_width,
	 error_message)
```
Sample values:
```
j_nu_py(): j_nu_py(230e9, 30, 1, 1.047, symphonyPy.MAXWELL_JUETTNER, symphonyPy.STOKES_I,
                   10, 2.5, 1, 1000, 1e10, 3.5, 10)
```
* Note: All parameters with units are in CGS.
* Note: In `C`, the keys `symphonyPy.MAXWELL_JUETTNER` and `symphonyPy.STOKES_I` are members of a struct called `params`.  They can be used with: `params->MAXWELL_JUETTNER` and `params->STOKES_I`.

### TODO
1. Add an anisotropic DF
2. Put in warnings for frequencies outside the intended frequency regime (much greater than the plasma, relativistic cyclotron frequencies)
3. Make cmake automatically add kernel_samples_datafiles/ directory to susceptibilityPy.pyx
