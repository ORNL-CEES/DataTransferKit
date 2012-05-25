//---------------------------------------------------------------------------//
/*!
 * \file DTK_CoreTypes.hpp
 * \author Stuart R. Slattery
 * \brief Types for the core subpackage
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CORETYPES_HPP
#define DTK_CORETYPES_HPP

namespace DataTransferKit
{

// Element topology enumerations.
enum DTK_ElementTopology {
    DTK_ElementTopology_MIN = 0,
    DTK_POINT = DTK_ElementTopology_MIN,
    DTK_LINE_SEGMENT,
    DTK_TRIANGLE,
    DTK_QUADRILATERAL,
    DTK_TETRAHEDRON,
    DTK_HEXAHEDRON,
    DTK_WEDGE,
    DTK_ElementTopology_MAX = DTK_WEDGE
};

} // end namespace DataTransferKit

#endif // end DTK_CORETYPES_HPP

//---------------------------------------------------------------------------//
// end DTK_CoreTypes.hpp
//---------------------------------------------------------------------------//

