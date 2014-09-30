##---------------------------------------------------------------------------##
## DataTransferKit/CMakeLists.txt
## Stuart R. Slattery
##---------------------------------------------------------------------------##
## TriBITS package file for the DataTransferKit package
##---------------------------------------------------------------------------##

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
