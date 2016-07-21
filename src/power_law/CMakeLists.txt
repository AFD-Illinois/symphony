add_library(power_law 
            ../distribution_function_common_routines.c 
            ../distribution_function_common_routines.h
            ../params.c
            ../params.h
            power_law.c 
            power_law.h 
            power_law_fits.c)
target_link_libraries(power_law ${GSL_LIBRARIES})
