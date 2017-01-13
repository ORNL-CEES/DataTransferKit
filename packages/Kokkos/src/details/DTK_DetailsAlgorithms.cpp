#include <details/DTK_DetailsAlgorithms.hpp>

#include <cmath>

namespace DataTransferKit
{
namespace Details
{

double distance( Point const &a, Point const &b )
{
    double distance_squared = 0.0;
    for ( int d = 0; d < 3; ++d )
    {
        double tmp = b[d] - a[d];
        distance_squared += tmp * tmp;
    }
    return std::sqrt( distance_squared );
}

double distance( Point const &point, Box const &box )
{
    Point projected_point;
    for ( int d = 0; d < 3; ++d )
    {
        if ( point[d] < box[2 * d + 0] )
            projected_point[d] = box[2 * d + 0];
        else if ( point[d] > box[2 * d + 1] )
            projected_point[d] = box[2 * d + 1];
        else
            projected_point[d] = point[d];
    }
    return distance( point, projected_point );
}

void expand( Box &box, Point const &point )
{
    for ( int d = 0; d < 3; ++d )
    {
        if ( point[d] < box[2 * d + 0] )
            box[2 * d + 0] = point[d];
        if ( point[d] > box[2 * d + 1] )
            box[2 * d + 1] = point[d];
    }
}

void expand( Box &box, Box const &other )
{
    for ( int d = 0; d < 3; ++d )
    {
        if ( box[2 * d + 0] > other[2 * d + 0] )
            box[2 * d + 0] = other[2 * d + 0];
        if ( box[2 * d + 1] < other[2 * d + 1] )
            box[2 * d + 1] = other[2 * d + 1];
    }
}

bool overlaps( Box const &box, Box const &other )
{
    for ( int d = 0; d < 3; ++d )
        if ( box[2 * d + 0] > other[2 * d + 1] ||
             box[2 * d + 1] < other[2 * d + 0] )
            return false;
    return true;
}

void centroid( Box const &box, Point &c )
{
    for ( int d = 0; d < 3; ++d )
        c[d] = 0.5 * ( box[2 * d + 0] + box[2 * d + 1] );
}

} // end namespace Details
} // end namespace DataTransferKit
