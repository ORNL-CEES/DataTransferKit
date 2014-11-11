SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Nanoflann             packages/Nanoflann              SS  REQUIRED
  Triangle              packages/Triangle               SS  OPTIONAL
  Utils                 packages/Utils                  SS  REQUIRED
  Interface             packages/Interface              SS  REQUIRED
  BasicGeometryAdapters packages/Adapters/BasicGeometry SS  OPTIONAL
  Operators             packages/Operators              SS  REQUIRED
  PointCloud            packages/PointCloud             SS  OPTIONAL
  STKMeshAdapters       packages/Adapters/STKMesh       SS  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
