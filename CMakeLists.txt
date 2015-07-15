##---------------------------------------------------------------------------##
## DataTransferKit/CMakeLists.txt
## Stuart R. Slattery
##---------------------------------------------------------------------------##
## TriBITS package file for the DataTransferKit package
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_DECL(DataTransferKit)
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_ETI_INSTANTIATIONS(TpetraCore
  "S={double} LO={int} GO={unsigned long} N={${KokkosClassic_DefaultNode}}"
  )

##---------------------------------------------------------------------------##
## Set the path to the DTKData directory.
##---------------------------------------------------------------------------##

SET( DTK_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DTKData" CACHE INTERNAL "" )
MESSAGE("DTK_DATA_DIR = '${DTK_DATA_DIR}'")

##---------------------------------------------------------------------------##
## Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PROCESS_SUBPACKAGES()
TRIBITS_PACKAGE_DEF()
TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
## Print the revision number to stdout    
##---------------------------------------------------------------------------##

FIND_PACKAGE(Git)
IF(GIT_FOUND)
    EXECUTE_PROCESS(
        COMMAND           ${GIT_EXECUTABLE} log --pretty=format:%H -n 1
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE   DTK_REVISION
        )
    MESSAGE("DTK Revision = '${DTK_REVISION}'")
ENDIF()

IF(Trilinos_ENABLE_DataTransferKitMoabAdapters)
    SET(MOAB_STAMP_FILE ${MOAB_LIBRARY_DIRS}/moab-revision)
    EXECUTE_PROCESS(COMMAND touch ${MOAB_STAMP_FILE})
    FILE(READ ${MOAB_STAMP_FILE} MOAB_REVISION)
    STRING(REGEX REPLACE "\n" "" MOAB_REVISION "${MOAB_REVISION}")
    MESSAGE("MOAB Revision = '${MOAB_REVISION}'")
ENDIF()

##---------------------------------------------------------------------------##
## end of DataTransferKit/CMakeLists.txt
##---------------------------------------------------------------------------##
