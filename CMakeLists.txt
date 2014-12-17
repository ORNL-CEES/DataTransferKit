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

TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra
  "S={double} LO={int} GO={unsigned long} N={${KokkosClassic_DefaultNode}}"
  "S={unsigned long} LO={int} GO={unsigned long} N={${KokkosClassic_DefaultNode}}"
  "S={int} LO={int} GO={unsigned long} N={${KokkosClassic_DefaultNode}}"
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
## end of DataTransferKit/CMakeLists.txt
##---------------------------------------------------------------------------##
