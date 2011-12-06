INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LG_Indexer_test
  SOURCES tstLG_Indexer.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Message_Buffer_test
#   SOURCES tstMessage_Buffer.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Interfaces_test
#   SOURCES tstInterfaces.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Transfer_Data_Field_test
#   SOURCES tstTransfer_Data_Field.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Transfer_Map_test
#   SOURCES tstTransfer_Map.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 4
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Mapper_test
#   SOURCES tstMapper.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 4
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Messenger_test
#   SOURCES tstMessenger.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 4
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Data_Transfer_Manager_test
#   SOURCES tstData_Transfer_Manager.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 4
#   STANDARD_PASS_OUTPUT
#   )
