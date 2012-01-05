##---------------------------------------------------------------------------##
## coupler/CMakeLists.txt
## Stuart R. Slattery
## Tuesday Dec 6 2011
##---------------------------------------------------------------------------##
## TriBITS package file for the Coupler package
##---------------------------------------------------------------------------##

INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(CombinedOption)

##---------------------------------------------------------------------------##
## A) Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_DECL(Coupler)

##---------------------------------------------------------------------------##
## B) Set up package-specific options
##---------------------------------------------------------------------------##
TRIBITS_ADD_DEBUG_OPTION()

##---------------------------------------------------------------------------##
## C) Add the libraries, tests, and examples
##---------------------------------------------------------------------------##

SET(${PACKAGE_NAME}_INCLUDE_DIRS "")
SET(${PACKAGE_NAME}_LIBRARY_DIRS "")
SET(${PACKAGE_NAME}_LIBRARIES "")

ADD_SUBDIRECTORY(src)
TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

##---------------------------------------------------------------------------##
## F) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_DEF()

TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of coupler/CMakeLists.txt
##---------------------------------------------------------------------------##
