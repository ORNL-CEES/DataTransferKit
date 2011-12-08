INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   LG_Indexer_test
#   SOURCES tstLG_Indexer.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Message_Buffer_test
#   SOURCES tstMessage_Buffer.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Transfer_Map_test
#   SOURCES tstTransfer_Map.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Interfaces_test
  SOURCES tstInterfaces.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Transfer_Data_Field_test
  SOURCES tstTransfer_Data_Field.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Mapper_test
  SOURCES tstMapper.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Messenger_test
#   SOURCES tstMessenger.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   Data_Transfer_Manager_test
#   SOURCES tstData_Transfer_Manager.cc ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   )
