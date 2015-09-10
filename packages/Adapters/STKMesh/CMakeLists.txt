##---------------------------------------------------------------------------##
## DataTransferKit/CMakeLists.txt
## Stuart R. Slattery
##---------------------------------------------------------------------------##
## TriBITS package file for the DataTransferKit package
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE(STKMeshAdapters)

##---------------------------------------------------------------------------##
## Add library, test, and examples.
##---------------------------------------------------------------------------##

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

IF( NOT ${PACKAGE_NAME}_FIND_DTK_DATA_DIR STREQUAL ${PACKAGE_NAME}_FIND_DTK_DATA_DIR-NOTFOUND)
    TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
ENDIF()

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
## end of DataTransferKit/CMakeLists.txt
##---------------------------------------------------------------------------##
