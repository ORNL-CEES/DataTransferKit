#ifndef DTK_DETAILS_ALGORITHMS_HPP
#define DTK_DETAILS_ALGORITHMS_HPP

#include <DTK_LinearBVH.hpp> // AABB

namespace DataTransferKit
{
namespace Details
{
using Point = Kokkos::Array<double, 3>;
using Box = AABB;

// distance point-point
double distance( Point const &a, Point const &b );
// distance point-box
double distance( Point const &point, Box const &box );
// expand an axis-aligned bounding box to include a point
void expand( Box &box, Point const &point );
// expand an axis-aligned bounding box to include another box
void expand( Box &box, Box const &other );
// check if two axis-aligned bounding boxes overlap
bool overlaps( Box const &box, Box const &other );
// calculate the centroid of a box
void centroid( Box const &box, Point &c );

} // end namespace Details
} // end namespace DataTransferKit

#endif
