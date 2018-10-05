IF (Trilinos_ENABLE_Fortran)
    # Enable unlimited-length lines and C-like preprocessing
    SET(CMAKE_Fortran_FLAGS "-cpp ${CMAKE_Fortran_FLAGS}")
ENDIF()


ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
