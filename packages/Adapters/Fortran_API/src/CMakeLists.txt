TRIBITS_ADD_LIBRARY(
  dtk_fortran_api
  HEADERS ""
  SOURCES DTK_Fortran_API.F90
  )
TARGET_LINK_LIBRARIES(dtk_fortran_api dtk_c_api)
