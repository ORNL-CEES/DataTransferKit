INCLUDE(TribitsTplDeclareLibraries)

TRIBITS_TPL_DECLARE_LIBRARIES( MOAB
  REQUIRED_HEADERS  iBase.h iMesh.h iMesh_extensions.h
  REQUIRED_LIBS_NAMES MOAB iMesh
  )
