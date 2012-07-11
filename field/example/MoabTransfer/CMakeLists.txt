INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MoabTransfer_1
  SOURCES moab_transfer_1.cpp PeaksEvaluator MoabMesh ArrayField
  COMM serial mpi
  DEPLIBS teuchos Utils Mesh Field
  )
