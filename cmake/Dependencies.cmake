TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Utils                 packages/Utils                  SS  REQUIRED
    Interface             packages/Interface              SS  REQUIRED
    BasicGeometryAdapters packages/Adapters/BasicGeometry SS  OPTIONAL
    Operators             packages/Operators              SS  REQUIRED
    IntrepidAdapters      packages/Adapters/Intrepid      SS  OPTIONAL
    STKMeshAdapters       packages/Adapters/STKMesh       SS  OPTIONAL
    MoabAdapters          packages/Adapters/Moab          SS  OPTIONAL
    LibmeshAdapters       packages/Adapters/Libmesh       SS  OPTIONAL
    ClassicDTKAdapters    packages/Adapters/ClassicDTK    SS  OPTIONAL
    PlainOldData          packages/Adapters/POD           SS  OPTIONAL
  )
