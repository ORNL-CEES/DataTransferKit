TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Utils                 packages/Utils                      ST  REQUIRED
    Interface             packages/Interface                  ST  OPTIONAL
    Search                packages/Search                     ST  REQUIRED
    Discretization        packages/Discretization             ST  OPTIONAL
    Meshfree              packages/Meshfree                   ST  OPTIONAL
    MapFactory            packages/MapFactory                 ST  OPTIONAL
    HybridTransport       packages/Benchmarks/HybridTransport ST  OPTIONAL
  )
