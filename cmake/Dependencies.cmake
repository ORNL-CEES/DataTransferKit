TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Utils                 packages/Utils                  SS  REQUIRED
    Kokkos                packages/Kokkos                 SS  REQUIRED
    Interface             packages/Interface              SS  REQUIRED
    BasicGeometryAdapters packages/Adapters/BasicGeometry SS  OPTIONAL
    IntrepidAdapters      packages/Adapters/Intrepid      SS  OPTIONAL
    Operators             packages/Operators              SS  REQUIRED
    C_API                 packages/Adapters/C_API         SS  OPTIONAL
    Fortran_API           packages/Adapters/Fortran_API   SS  OPTIONAL
    STKMeshAdapters       packages/Adapters/STKMesh       SS  OPTIONAL
    MoabAdapters          packages/Adapters/Moab          SS  OPTIONAL
    LibmeshAdapters       packages/Adapters/Libmesh       SS  OPTIONAL
    ClassicDTKAdapters    packages/Adapters/ClassicDTK    SS  OPTIONAL
  )
