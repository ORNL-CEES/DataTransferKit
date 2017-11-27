# While DTK does not depend on Shards Intrepid2 does so we need to add Shards as
# a dependency
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
