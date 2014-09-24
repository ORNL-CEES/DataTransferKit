INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_BoundingBox.hpp
  DTK_Box.hpp
  DTK_BuhmannBasis.hpp
  DTK_BuhmannBasis_impl.hpp
  DTK_CellTopologyFactory.hpp
  DTK_CenterDistributor.hpp
  DTK_CenterDistributor_impl.hpp
  DTK_CloudDomain.hpp
  DTK_CloudDomain_impl.hpp
  DTK_CommIndexer.hpp
  DTK_CommTools.hpp
  DTK_Cylinder.hpp
  DTK_DBC.hpp
  DTK_ElementMeasure.hpp
  DTK_ElementTree.hpp
  DTK_EuclideanDistance.hpp
  DTK_EuclideanDistance_impl.hpp
  DTK_FieldContainer.hpp
  DTK_FieldEvaluator.hpp
  DTK_FieldIntegrator.hpp
  DTK_FieldManager.hpp
  DTK_FieldManager_def.hpp
  DTK_FieldTools.hpp
  DTK_FieldTools_def.hpp
  DTK_FieldTraits.hpp
  DTK_GeometryManager.hpp
  DTK_GeometryManager_def.hpp
  DTK_GeometryRCB.hpp
  DTK_GeometryRCB_def.hpp
  DTK_GeometryRendezvous.hpp
  DTK_GeometryRendezvous_def.hpp
  DTK_GeometryTraits.hpp
  DTK_IntegralAssemblyMap.hpp
  DTK_IntegralAssemblyMap_def.hpp
  DTK_IntrepidCell.hpp
  DTK_IntrepidCell_impl.hpp
  DTK_IntrepidSideCell.hpp
  DTK_IntrepidSideCell_impl.hpp
  DTK_LocalMLSProblem.hpp
  DTK_LocalMLSProblem_impl.hpp
  DTK_MeshBlock.hpp
  DTK_MeshContainer.hpp
  DTK_MeshFreeInterpolator.hpp
  DTK_MeshFreeInterpolatorFactory.hpp
  DTK_MeshFreeInterpolatorFactory_impl.hpp
  DTK_MeshManager.hpp
  DTK_MeshTools.hpp
  DTK_MeshTraits.hpp
  DTK_MeshTraitsFieldAdapter.hpp
  DTK_MeshTypes.hpp
  DTK_MovingLeastSquare.hpp
  DTK_MovingLeastSquare_impl.hpp
  DTK_Partitioner.hpp
  DTK_PartitionerFactory.hpp
  DTK_PolynomialMatrix.hpp
  DTK_PolynomialMatrix_impl.hpp
  DTK_RadialBasisPolicy.hpp
  DTK_RCB.hpp
  DTK_Rendezvous.hpp
  DTK_Rendezvous_def.hpp
  DTK_SearchTreeFactory.hpp
  DTK_SerialPartitioner.hpp
  DTK_SharedDomainMap.hpp
  DTK_SharedDomainMap_def.hpp
  DTK_SplineInterpolationPairing.hpp
  DTK_SplineInterpolationPairing_impl.hpp
  DTK_SplineInterpolator.hpp
  DTK_SplineInterpolator_impl.hpp
  DTK_SplineOperatorA.hpp
  DTK_SplineOperatorA_impl.hpp
  DTK_SplineOperatorC.hpp
  DTK_SplineOperatorC_impl.hpp
  DTK_StaticSearchTree.hpp
  DTK_StaticSearchTree_impl.hpp
  DTK_TopologyTools.hpp
  DTK_TopologyTools_def.hpp
  DTK_VolumeSourceMap.hpp
  DTK_VolumeSourceMap_def.hpp
  DTK_WendlandBasis.hpp
  DTK_WendlandBasis_impl.hpp
  DTK_WuBasis.hpp
  DTK_WuBasis_impl.hpp
  ) 

APPEND_SET(SOURCES
  DTK_BoundingBox.cpp
  DTK_Box.cpp
  DTK_CellTopologyFactory.cpp
  DTK_CommIndexer.cpp
  DTK_CommTools.cpp
  DTK_Cylinder.cpp
  DTK_DBC.cpp
  DTK_ElementTree.cpp
  DTK_MeshManager.cpp
  DTK_MeshTools.cpp
  DTK_RCB.cpp
  DTK_SearchTreeFactory.cpp
  DTK_SerialPartitioner.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  datatransferkit
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
