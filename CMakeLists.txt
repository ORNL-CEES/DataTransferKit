##---------------------------------------------------------------------------##
## DataTransferKit/CMakeLists.txt
## Stuart R. Slattery
##---------------------------------------------------------------------------##
## TriBITS package file for the DataTransferKit package
##---------------------------------------------------------------------------##

INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(CombinedOption)

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(DataTransferKit)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##

# Design-by-Contract
TRIBITS_ADD_OPTION_AND_DEFINE(
	DataTransferKit_ENABLE_DBC
	HAVE_DTK_DBC
	"Enable Design-by-Contract checks. WARNING: use for debug checking but disable for production runs as it incurs a significant overhead."
	${${PROJECT_NAME}_ENABLE_DEBUG}
)

# If Zoltan and MPI must BOTH be enabled to function in parallel. Therefore 
# here we turn off MPI support for DataTransferKit explicitly if both are
# not enabled.
IF( ${PACKAGE_NAME}_ENABLE_MPI AND ${PROJECT_NAME}_ENABLE_Zoltan )
  GLOBAL_SET( HAVE_DTK_MPI ON )
ELSE()
  IF (${PACKAGE_NAME}_ENABLE_MPI)
    MESSAGE("WARNING: ${PACKAGE_NAME}_ENABLE_MPI=TRUE but ${PROJECT_NAME}_ENABLE_Zoltan=FALSE,"
     " DataTransferKit can't provide parallel functionality!")
  ENDIF()
  GLOBAL_SET( HAVE_DTK_MPI OFF )
ENDIF()

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

APPEND_SET(Tpetra_ETI_GORDS "unsigned int")
GLOBAL_SET(Tpetra_ETI_GORDS ${Tpetra_ETI_GORDS})
APPEND_SET(Tpetra_ETI_LORDS "unsigned int")
GLOBAL_SET(Tpetra_ETI_LORDS ${Tpetra_ETI_LORDS})

##---------------------------------------------------------------------------##
## C) Add library, test, and examples.
##---------------------------------------------------------------------------##

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of DataTransferKit/CMakeLists.txt
##---------------------------------------------------------------------------##
