GLOBAL_SET(BoostOrg_INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
GLOBAL_SET(BoostOrg_LIBRARY_DIRS "${Boost_LIBRARY_DIRS}")

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES(
  BoostOrg
  REQUIRED_HEADERS boost/property_tree/json_parser.hpp
                   boost/property_tree/ptree.hpp
                   boost/math/tools/polynomial.hpp
                   boost/math/tools/rational.hpp
                   boost/geometry.hpp
                   boost/range.hpp
                   boost/program_options.hpp
                   boost/test/unit_test.hpp
  )

# Use CMake FindBoost module to check version is sufficient
SET(BOOST_INCLUDEDIR ${TPL_BoostOrg_INCLUDE_DIRS})
MESSAGE(STATUS "BOOST_INCLUDEDIR: ${BOOST_INCLUDEDIR}")
SET(Boost_NO_SYSTEM_PATHS ON)
FIND_PACKAGE(Boost 1.61.0 REQUIRED)
