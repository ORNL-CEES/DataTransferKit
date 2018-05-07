TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Utils                 packages/Utils                      ST  REQUIRED
    Interface             packages/Interface                  ST  REQUIRED
    Search                packages/Search                     ST  REQUIRED
    Discretization        packages/Discretization             ST  REQUIRED
    Meshfree              packages/Meshfree                   ST  REQUIRED
    HybridTransport       packages/Benchmarks/HybridTransport ST  OPTIONAL
  )
