##---------------------------------------------------------------------------##
## DataTransferKit/CMakeLists.txt
## Stuart R. Slattery
##---------------------------------------------------------------------------##
## TriBITS package file for the DataTransferKit package
##---------------------------------------------------------------------------##

INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_DECL(DataTransferKit)
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

##---------------------------------------------------------------------------##
## Set the path to the DTKData directory.
##---------------------------------------------------------------------------##

SET( DTKDATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DTKData" CACHE PATH "" )

##---------------------------------------------------------------------------##
## Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PROCESS_SUBPACKAGES()
TRIBITS_PACKAGE_DEF()
TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
## end of DataTransferKit/CMakeLists.txt
##---------------------------------------------------------------------------##
