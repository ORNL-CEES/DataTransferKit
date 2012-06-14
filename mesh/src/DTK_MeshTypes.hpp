//---------------------------------------------------------------------------//
/*!
 * \file DTK_MeshTypes.hpp
 * \author Stuart R. Slattery
 * \brief Enumerated types for the mesh subpackage.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTYPES_HPP
#define DTK_MESHTYPES_HPP

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
    DTK_PYRAMID,
    DTK_ElementTopology_MAX = DTK_PYRAMID
};

} // end namespace DataTransferKit

#endif // end DTK_MESHTYPES_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshTypes.hpp
//---------------------------------------------------------------------------//

