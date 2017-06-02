##---------------------------------------------------------------------------##
## DataTransferKit/CMakeLists.txt
## Stuart R. Slattery
##---------------------------------------------------------------------------##
## TriBITS package file for the DataTransferKit package
##---------------------------------------------------------------------------##

INCLUDE(CombinedOption)
INCLUDE(TribitsETISupport)

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE(Utils)

# Design-by-Contract
TRIBITS_ADD_OPTION_AND_DEFINE(
        DataTransferKit_ENABLE_DBC
        HAVE_DTK_DBC
        "Enable Design-by-Contract checks. WARNING: use for debug checking but disable for production runs as it incurs a significant overhead."
        ${${PROJECT_NAME}_ENABLE_DEBUG}
)

# ETI
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_ETI_SUPPORT()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)

TRIBITS_ADD_OPTION_AND_DEFINE(
  DataTransferKit_ENABLE_EXPLICIT_INSTANTIATION
  HAVE_DATATRANSFERKIT_EXPLICIT_INSTANTIATION
  "Enable explicit template instantiation (ETI) in DTK"
  ${${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
  )

##---------------------------------------------------------------------------##
## Add library, test, and examples.
##---------------------------------------------------------------------------##

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
##                      end of DataTransferKit/CMakeLists.txt
##---------------------------------------------------------------------------##
