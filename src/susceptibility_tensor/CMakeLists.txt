add_library(susceptibility_tensor susceptibility_tensor.h mueller_matrix.c chi_ij.c integrands_suscept.c integrator.c ../params.c ../params.h)
target_link_libraries(susceptibility_tensor ${GSL_LIBRARIES})

cython_add_module(susceptibilityPy susceptibilityPy.pyx)

target_link_libraries(susceptibilityPy susceptibility_tensor)
