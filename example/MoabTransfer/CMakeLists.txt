INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MoabTransfer_1
  SOURCES moab_transfer_1.cpp PeaksEvaluator MoabMesh ArrayField
  COMM serial mpi
  DEPLIBS datatransferkit
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  MoabTransferINP
  SOURCE_FILES tri_peaks.vtk quad_mesh.vtk
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS MoabTransfer_1
)
