SET(${PACKAGE_NAME}_Trilinos_REQUIRED_COMPONENTS Belos Kokkos Intrepid2 Stratimikos Teuchos Thyra Tpetra)
SET(${PACKAGE_NAME}_Trilinos_OPTIONAL_COMPONENTS "")

IF (${PACKAGE_NAME}_TRILINOS_TPL)
  TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
    LIB_REQUIRED_TPLS
      MPI
      Trilinos
      BoostOrg
    LIB_OPTIONAL_TPLS
      BoostOrg
      Netcdf
    TEST_OPTIONAL_TPLS
      # BoostOrg is listed twice to have -isystem added to the include directories for the TPL when compiling the tests.
      # This is a known limitation of TriBITS (c.f. https://tribits.org/doc/TribitsDevelopersGuide.html#project-name-tpl-system-include-dirs)
      BoostOrg
    )
ELSE()
  TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
    LIB_REQUIRED_PACKAGES
      ${${PACKAGE_NAME}_Trilinos_REQUIRED_COMPONENTS}
    LIB_OPTIONAL_PACKAGES
      ${${PACKAGE_NAME}_Trilinos_OPTIONAL_COMPONENTS}
    LIB_REQUIRED_TPLS
      MPI
    LIB_OPTIONAL_TPLS
      BoostOrg
      Netcdf
    TEST_OPTIONAL_TPLS
      # BoostOrg is listed twice to have -isystem added to the include directories for the TPL when compiling the tests.
      # This is a known limitation of TriBITS (c.f. https://tribits.org/doc/TribitsDevelopersGuide.html#project-name-tpl-system-include-dirs)
      BoostOrg
    )
ENDIF()
