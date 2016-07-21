add_library(maxwell_juettner
            ../distribution_function_common_routines.c 
            ../distribution_function_common_routines.h
            ../params.c
            ../params.h 
            maxwell_juettner.c 
            maxwell_juettner.h
            maxwell_juettner_fits.c)
target_link_libraries(maxwell_juettner ${GSL_LIBRARIES})
