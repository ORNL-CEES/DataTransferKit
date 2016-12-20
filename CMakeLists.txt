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

##---------------------------------------------------------------------------##
## Explicit template instantiation
##---------------------------------------------------------------------------##

# Setup the same ETI as Tpetra but do not set global ordinal types as
# additional scalar types.

# check tpetra variables are defined.
ASSERT_DEFINED(TpetraCore_ETI_SCALARS)
ASSERT_DEFINED(TpetraCore_ETI_LORDS)
ASSERT_DEFINED(TpetraCore_ETI_GORDS)
ASSERT_DEFINED(TpetraCore_ETI_NODES)

# set the eti types
SET (DataTransferKit_ETI_SCALARS_TEMP ) # start with the empty list
FOREACH (scalarName IN LISTS TpetraCore_ETI_SCALARS)
  LIST (FIND TpetraCore_ETI_GORDS ${scalarName} scalarNameIndex)
  IF (-1 EQUAL scalarNameIndex) # not in the list of integer types
    LIST (APPEND DataTransferKit_ETI_SCALARS_TEMP ${scalarName})
  ENDIF ()
ENDFOREACH ()

SET(DataTransferKit_ETI_SCALARS ${DataTransferKit_ETI_SCALARS_TEMP})
SET(DataTransferKit_ETI_LORDS   ${TpetraCore_ETI_LORDS})
SET(DataTransferKit_ETI_GORDS   ${TpetraCore_ETI_GORDS})
SET(DataTransferKit_ETI_NODES   ${TpetraCore_ETI_NODES})

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_ETI_SUPPORT()

# make the eti variables globally accessible
GLOBAL_SET(DataTransferKit_ETI_SCALARS ${DataTransferKit_ETI_SCALARS})
GLOBAL_SET(DataTransferKit_ETI_LORDS   ${DataTransferKit_ETI_LORDS})
GLOBAL_SET(DataTransferKit_ETI_GORDS   ${DataTransferKit_ETI_GORDS})
GLOBAL_SET(DataTransferKit_ETI_NODES   ${DataTransferKit_ETI_NODES})

##---------------------------------------------------------------------------##
## Set the path to the DTKData directory.
##---------------------------------------------------------------------------##

# Search for one of the mesh files in the DTKData directory
find_file(DTK_DATA_MESH_FILE cube_mesh.exo PATHS
    ${${PACKAGE_NAME}_SOURCE_DIR}/data/exodus
    ${${PACKAGE_NAME}_SOURCE_DIR}/DTKData/exodus
    NO_DEFAULT_PATH
)
if(DTK_DATA_MESH_FILE)
    # Get the location of the mesh file
    get_filename_component(DTK_DATA_DIR ${DTK_DATA_MESH_FILE} DIRECTORY)
    # Strip "/exodus" from the path
    string(REGEX REPLACE "/exodus" "" DTK_DATA_DIR ${DTK_DATA_DIR})
    message("-- DTKData found at ${DTK_DATA_DIR}")
else()
    set(DTK_DATA_DIR DTK_DATA_DIR-NOTFOUND)
    message("-- DTKData NOT found")
endif()

##---------------------------------------------------------------------------##
## Disable Fortran_API if Trilinos_ENABLE_Fortran is FALSE
##---------------------------------------------------------------------------##
IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    IF (${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}Fortran_API})
        MESSAGE("-- " "NOTE: Setting ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}Fortran_API=OFF because"
                "${PROJECT_NAME}_ENABLE_Fortran = '${${PROJECT_NAME}_ENABLE_Fortran}'")
        SET(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}Fortran_API OFF)
    ENDIF()
ENDIF()

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

##---------------------------------------------------------------------------##
## Check C++ code style using clang-format version 3.9
##---------------------------------------------------------------------------##

IF(DataTransferKit_ENABLE_ClangFormat)
    INCLUDE(CodeFormat)
ENDIF()

##---------------------------------------------------------------------------##
## Check that there are no trailing whitespaces or tabs
##---------------------------------------------------------------------------##

INCLUDE(NoTrailingWhitespacesOrTabs)

##---------------------------------------------------------------------------##
## Build the documentation published on http://datatransferkit.readthedocs.org
##---------------------------------------------------------------------------##

IF(DataTransferKit_ENABLE_ReadTheDocs)
    ADD_SUBDIRECTORY(docs)
ENDIF()

##---------------------------------------------------------------------------##
## Generate the Doxygen documentation
##---------------------------------------------------------------------------##

IF(DataTransferKit_ENABLE_Doxygen)
    ADD_SUBDIRECTORY(docs/doxygen)
ENDIF()

##---------------------------------------------------------------------------##
## Generate a .ycm_extra_conf.py file for use with YouCompleteMe
##---------------------------------------------------------------------------##

IF(DataTransferKit_ENABLE_YouCompleteMe)
    INCLUDE(CodeCompletion)
ENDIF()

##---------------------------------------------------------------------------##
## end of DataTransferKit/CMakeLists.txt
##---------------------------------------------------------------------------##
