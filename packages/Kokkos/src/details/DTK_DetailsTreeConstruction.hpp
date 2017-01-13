#ifndef DTK_DETAILS_TREE_CONSTRUCTION_HPP
#define DTK_DETAILS_TREE_CONSTRUCTION_HPP

#include <DTK_LinearBVH.hpp>
#include <Kokkos_Pair.hpp>

namespace DataTransferKit
{
namespace Details
{
// utilities for tree construction
unsigned int expandBits( unsigned int v );
unsigned int morton3D( double x, double y, double z );
int countLeadingZeros( unsigned int k );
int commonPrefix( unsigned int const *k, int n, int i, int j );
int findSplit( unsigned int *sorted_morton_codes, int first, int last );
Kokkos::pair<int, int> determineRange( unsigned int *sorted_morton_codes, int n,
                                       int i );
// COMMENT: most of these could/should be protected function in BVH to avoid
// passing all this data around
void calculateBoundingBoxOfTheScene( AABB const *bounding_boxes, int n,
                                     AABB &scene_bounding_box );
void assignMortonCodes( AABB const *bounding_boxes, unsigned int *morton_codes,
                        int n, AABB const &scene_bounding_box );
void sortObjects( unsigned int *morton_codes, int *object_ids, int n );
Node *generateHierarchy( unsigned int *sorted_morton_codes, int n,
                         Node *leaf_nodes, Node *internal_nodes );
void calculateBoundingBoxes( Node const *leaf_nodes, Node *internal_nodes,
                             int n );

} // end namespace Details
} // end namespace DataTransferKit

#endif
