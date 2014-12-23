
#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_BuhmannBasis.hpp
  DTK_BuhmannBasis_impl.hpp
  DTK_CenterDistributor.hpp
  DTK_CenterDistributor_impl.hpp
  DTK_CloudDomain.hpp
  DTK_CloudDomain_impl.hpp
  DTK_EuclideanDistance.hpp
  DTK_EuclideanDistance_impl.hpp
  DTK_LocalMLSProblem.hpp
  DTK_LocalMLSProblem_impl.hpp
  DTK_MovingLeastSquareReconstructionOperator.hpp
  DTK_MovingLeastSquareReconstructionOperator_impl.hpp
  DTK_PointCloudDummy.hpp
  DTK_PolynomialMatrix.hpp
  DTK_PolynomialMatrix_impl.hpp
  DTK_RadialBasisPolicy.hpp
  DTK_SplineCoefficientMatrix.hpp
  DTK_SplineCoefficientMatrix_impl.hpp
  DTK_SplineEvaluationMatrix.hpp
  DTK_SplineEvaluationMatrix_impl.hpp
  DTK_SplineInterpolationPairing.hpp
  DTK_SplineInterpolationPairing_impl.hpp
  DTK_SplineInterpolationOperator.hpp
  DTK_SplineInterpolationOperator_impl.hpp
  DTK_SplineProlongationOperator.hpp
  DTK_SplineProlongationOperator_impl.hpp
  DTK_WendlandBasis.hpp
  DTK_WendlandBasis_impl.hpp
  DTK_WuBasis.hpp
  DTK_WuBasis_impl.hpp
  ) 

APPEND_SET(SOURCES
  DTK_PointCloudDummy.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_pointcloud
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
