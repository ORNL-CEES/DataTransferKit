#include <details/DTK_DetailsAlgorithms.hpp>

#include <cmath>

namespace DataTransferKit
{
namespace Details
{

void expand( BBox &box, Point const &point )
{
    for ( int d = 0; d < 3; ++d )
    {
        if ( point[d] < box[2 * d + 0] )
            box[2 * d + 0] = point[d];
        if ( point[d] > box[2 * d + 1] )
            box[2 * d + 1] = point[d];
    }
}

} // end namespace Details
} // end namespace DataTransferKit
