TODO: add __init__.py to build/ folder and build/susceptibility_tensor/ folder

#*symphony*

##Overview

*symphony* calculates synchrotron emissivities and absorptivities, polarized in the Stokes basis {I, Q, U, V}, for any arbitrary gyrotropic momentum space distribution function.  Three distribution functions are built in: a relativistic thermal (Maxwell-Juettner) distribution, a nonthermal power-law distribution, and a kappa distribution.

*symphony* is described in [Pandya et al., 2016 ApJ 822 34](http://dx.doi.org/10.3847/0004-637X/822/1/34). If you use this code in an academic context, cite that paper.

##Features

* `C` code to calculate synchrotron emissivities via the `j_nu()` function and absorptivities via the `alpha_nu()` function.
* `C` code to evaluate approximate fitting function values for the emissivity and absorptivity, via the `j_nu_fit()` and `alpha_nu_fit()` functions, respectively.
* CMake configure system, which helps during the build process to find all necessary libraries and files.
* `Python` interface for `j_nu()`, `alpha_nu()`, `j_nu_fit()`, and `alpha_nu_fit()`.
  * This combines the speed of `C` when evaluating emissivities and absorptivities with `Python`'s user-friendly syntax.  It also allows for interfacing with larger `Python` codes.

##How to use *symphony*

###To download and build *symphony*:
 1. Clone *symphony* from github.  Navigate into the *symphony* folder, and note that it contains a folder named "src/"; create a folder named "build" and navigate into it.
 2. Type "cmake" followed by the location of the "src/" folder.  Altogether, this line should look something like: "cmake /location/to/symphony/src". You can add the argument `-DCMAKE_INSTALL_PREFIX=/name/of/dir` to set the name of the directory to install to
 3. Type "make"
 4. Optionally, run `make install` to install the library and Python module onto your system.

###To use *symphony*'s `Python` interface:
 1. Navigate to the "build/" folder created in step 1., above.  Open `Python` in the command line or by writing a ".py" file.
 2. Import *symphony* by typing `import symphonyPy`.  
  * This allows one to call 4 functions: `j_nu_py()`, `alpha_nu_py()`, `j_nu_fit_py()`, and `alpha_nu_fit_py()`.  
  * The first two provide calculated values of the emissivity and absorptivity for the input parameters, and the latter two provide the corresponding approximate fitting formula results.
 3. The arguments of these functions can be found by accessing the associated docstrings.  This can be done in the `Python` command line using the following: 
```
import symphonyPy
symphonyPy.j_nu_py?
```

###Arguments for Emissivity and Absorptivity Functions
* Arguments for all emissivity and absorptivity functions in both `C` and `Python` take nearly the same arguments and output a double.  The only difference is that the C version has an `error_message` parameter for handling evaluation errors. The arguments are: 
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
