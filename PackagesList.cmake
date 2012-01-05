#
# See documentation in Trilinos preCopyrightTrilinos/ExtraExternalRepositories.cmake
#

INCLUDE(TribitsListHelpers)

SET( Coupler_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  Coupler         .     SS
  )

PACKAGE_DISABLE_ON_PLATFORMS(Coupler Windows)
