add_library(kappa
            ../distribution_function_common_routines.c 
            ../distribution_function_common_routines.h
            ../params.c
            ../params.h
            kappa.c 
            kappa.h 
            kappa_fits.c)
target_link_libraries(kappa ${GSL_LIBRARIES})
