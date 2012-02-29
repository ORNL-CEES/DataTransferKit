#
# See documentation in Trilinos preCopyrightTrilinos/ExtraExternalRepositories.cmake
#

INCLUDE(TribitsListHelpers)

SET( DataTransferKit_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  DataTransferKit         .     SS
  )

PACKAGE_DISABLE_ON_PLATFORMS(DataTransferKit Windows)
