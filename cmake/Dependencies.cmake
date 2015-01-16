TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Nanoflann             packages/Nanoflann              SS  REQUIRED
    Utils                 packages/Utils                  SS  REQUIRED
    Interface             packages/Interface              SS  REQUIRED
    BasicGeometryAdapters packages/Adapters/BasicGeometry SS  OPTIONAL
    Operators             packages/Operators              SS  REQUIRED
    IntrepidAdapters      packages/Adapters/Intrepid      SS  OPTIONAL
    STKMeshAdapters       packages/Adapters/STKMesh       SS  OPTIONAL
  )
