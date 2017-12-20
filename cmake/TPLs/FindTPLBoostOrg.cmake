GLOBAL_SET(BoostOrg_INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
GLOBAL_SET(BoostOrg_LIBRARY_DIRS "${Boost_LIBRARY_DIRS}")

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES(
  BoostOrg
  REQUIRED_HEADERS boost/version.hpp boost/mpl/at.hpp
  REQUIRED_LIBS_NAMES boost_mpi boost_serialization
  )
