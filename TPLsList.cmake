IF (${PACKAGE_NAME}_TRILINOS_TPL)
  SET(TRILINOS_TPL
    Trilinos "cmake/TPLs/" SS
   )
ELSE()
  SET(TRILINOS_TPL "")
ENDIF()

IF (${PACKAGE_NAME}_ARBORX_TPL)
SET(ARBORX_TPL
  ArborX "cmake/TPLs/" SS
 )
ELSE()
SET(ARBORX_TPL "")
ENDIF()

TRIBITS_REPOSITORY_DEFINE_TPLS(
  BoostOrg  "cmake/TPLs/"                                       SS
  MPI       "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"     SS
  Netcdf    "${${PROJECT_NAME}_TRIBITS_DIR}/common_tpls/"       SS
  ${TRILINOS_TPL}
  ${ARBORX_TPL}
)
