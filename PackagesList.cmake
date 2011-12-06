#
# See documentation in Trilinos preCopyrightTrilinos/ExtraExternalRepositories.cmake
#

INCLUDE(TribitsListHelpers)

SET( Coupler_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  coupler         ./     SS
  )

PACKAGE_DISABLE_ON_PLATFORMS(coupler Windows)
