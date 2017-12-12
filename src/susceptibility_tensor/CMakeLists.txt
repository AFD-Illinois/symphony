add_library(susceptibility_tensor susceptibility_tensor.h mueller_matrix.c chi_ij.c integrands.c integrator.c ../params.c ../params.h)

cython_add_module(susceptibilityPy susceptibilityPy.pyx)

target_link_libraries(susceptibilityPy susceptibility_tensor ${GSL_LIBRARIES} ${CBLAS_LIBRARIES})
