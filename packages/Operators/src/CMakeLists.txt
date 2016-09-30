SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Search)
APPEND_SET(HEADERS
  ${DIR}/DTK_CoarseGlobalSearch.hpp
  ${DIR}/DTK_CoarseLocalSearch.hpp
  ${DIR}/DTK_FineLocalSearch.hpp
  ${DIR}/DTK_ParallelSearch.hpp
  )

APPEND_SET(SOURCES
  ${DIR}/DTK_CoarseGlobalSearch.cpp
  ${DIR}/DTK_CoarseLocalSearch.cpp
  ${DIR}/DTK_FineLocalSearch.cpp
  ${DIR}/DTK_ParallelSearch.cpp
  )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/SharedDomain)
APPEND_SET(HEADERS
  ${DIR}/DTK_ConsistentInterpolationOperator.hpp
  ${DIR}/DTK_L2ProjectionOperator.hpp
  )

APPEND_SET(SOURCES
  ${DIR}/DTK_ConsistentInterpolationOperator.cpp
  ${DIR}/DTK_L2ProjectionOperator.cpp
  )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/PointCloud)
APPEND_SET(HEADERS
  ${DIR}/DTK_BuhmannBasis.hpp
  ${DIR}/DTK_BuhmannBasis_impl.hpp
  ${DIR}/DTK_CenterDistributor.hpp
  ${DIR}/DTK_CenterDistributor_impl.hpp
  ${DIR}/DTK_CloudDomain.hpp
  ${DIR}/DTK_CloudDomain_impl.hpp
  ${DIR}/DTK_EuclideanDistance.hpp
  ${DIR}/DTK_EuclideanDistance_impl.hpp
  ${DIR}/DTK_LocalMLSProblem.hpp
  ${DIR}/DTK_LocalMLSProblem_impl.hpp
  ${DIR}/DTK_MovingLeastSquareReconstructionOperator.hpp
  ${DIR}/DTK_MovingLeastSquareReconstructionOperator_impl.hpp
  ${DIR}/DTK_NodeToNodeOperator.hpp
  ${DIR}/DTK_NodeToNodeOperator_impl.hpp
  ${DIR}/DTK_PolynomialMatrix.hpp
  ${DIR}/DTK_RadialBasisPolicy.hpp
  ${DIR}/DTK_SplineCoefficientMatrix.hpp
  ${DIR}/DTK_SplineCoefficientMatrix_impl.hpp
  ${DIR}/DTK_SplineEvaluationMatrix.hpp
  ${DIR}/DTK_SplineEvaluationMatrix_impl.hpp
  ${DIR}/DTK_SplineInterpolationPairing.hpp
  ${DIR}/DTK_SplineInterpolationPairing_impl.hpp
  ${DIR}/DTK_SplineInterpolationOperator.hpp
  ${DIR}/DTK_SplineInterpolationOperator_impl.hpp
  ${DIR}/DTK_SplineProlongationOperator.hpp
  ${DIR}/DTK_WendlandBasis.hpp
  ${DIR}/DTK_WendlandBasis_impl.hpp
  ${DIR}/DTK_WuBasis.hpp
  ${DIR}/DTK_WuBasis_impl.hpp
  )

APPEND_SET(SOURCES
  ${DIR}/DTK_CloudDomain.cpp
  ${DIR}/DTK_LocalMLSProblem.cpp
  ${DIR}/DTK_MovingLeastSquareReconstructionOperator.cpp
  ${DIR}/DTK_NodeToNodeOperator.cpp
  ${DIR}/DTK_PolynomialMatrix.cpp
  ${DIR}/DTK_SplineCoefficientMatrix.cpp
  ${DIR}/DTK_SplineEvaluationMatrix.cpp
  ${DIR}/DTK_SplineInterpolationOperator.cpp
  ${DIR}/DTK_SplineInterpolationPairing.cpp
  ${DIR}/DTK_SplineProlongationOperator.cpp
  )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Factory)
APPEND_SET(HEADERS
  ${DIR}/DTK_MapOperatorFactory.hpp
  ${DIR}/DTK_PointCloudOperatorFactory.hpp
  )

APPEND_SET(SOURCES
  ${DIR}/DTK_MapOperatorFactory.cpp
  ${DIR}/DTK_PointCloudOperatorFactory.cpp
  )

TRIBITS_ADD_LIBRARY(
  dtk_operators
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
