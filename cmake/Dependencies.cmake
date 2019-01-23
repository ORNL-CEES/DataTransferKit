TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES
    Kokkos
    Teuchos
    Tpetra
    Intrepid2
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
