#ifndef DTK_LINEAR_BVH_DECL_HPP
#define DTK_LINEAR_BVH_DECL_HPP

#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_Array.hpp>
#include <Kokkos_View.hpp>

#include <DTK_Box.hpp>
#include <DTK_Node.hpp>
#include <details/DTK_DetailsAlgorithms.hpp>
#include <details/DTK_Predicate.hpp>

#include "DTK_ConfigDefs.hpp"

namespace DataTransferKit
{
/**
 * Bounding Volume Hierarchy.
 */
template <typename NO>
struct BVH
{
  public:
    using DeviceType = typename NO::device_type;

    BVH( Kokkos::View<BBox const *, DeviceType> bounding_boxes );

    Kokkos::View<Node *, DeviceType> leaf_nodes;
    Kokkos::View<Node *, DeviceType> internal_nodes;
    /**
     *  Array of indices that sort the boxes used to construct the hierarchy.
     *  The leaf nodes are ordered so we need these to identify objects that
     * meet
     *  a predicate.
     */
    Kokkos::View<int *, DeviceType> indices;
};

} // end namespace DataTransferKit

#endif
