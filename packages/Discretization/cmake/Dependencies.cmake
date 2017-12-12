# The dependence on Shards can be removed once we update Trilinos. Right now it
# is necessary due to a bug in Intrepid, which has been fixed in develop
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES
  DataTransferKitInterface
  DataTransferKitSearch
  DataTransferKitUtils
  Intrepid2
  Kokkos
  Shards
  Teuchos
  Tpetra
  TEST_REQUIRED_PACKAGES
  DataTransferKitInterface
  )
