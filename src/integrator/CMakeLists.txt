add_library(integrator integrate.c integrate.h integrands.c ${CMAKE_SOURCE_DIR}/bessel_mod.c)
target_link_libraries(integrator kappa maxwell_juettner power_law)
