//---------------------------------------------------------------------------//
/*!
 * \file DTK_CoreTypes.hpp
 * \author Stuart R. Slattery
 * \brief Enumerated types for the core subpackage.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CORETYPES_HPP
#define DTK_CORETYPES_HPP

namespace DataTransferKit
{

//! Element type enumerations.
enum DTK_ElementType {
    DTK_ElementType_MIN = 0,
    DTK_POINT = DTK_ElementType_MIN,
    DTK_EDGE,
    DTK_FACE,
    DTK_REGION,
    DTK_ElementType_MAX = DTK_REGION
};

//! Element topology enumerations.
enum DTK_ElementTopology {
    DTK_ElementTopology_MIN = 0,
    DTK_VERTEX = DTK_ElementTopology_MIN,
    DTK_LINE_SEGMENT,
    DTK_TRIANGLE,
    DTK_QUADRILATERAL,
    DTK_TETRAHEDRON,
    DTK_HEXAHEDRON,
    DTK_ElementTopology_MAX = DTK_HEXAHEDRON
};

} // end namespace DataTransferKit

#endif // end DTK_CORETYPES_HPP

//---------------------------------------------------------------------------//
// end DTK_CoreTypes.hpp
//---------------------------------------------------------------------------//

