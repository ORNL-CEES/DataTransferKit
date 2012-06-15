//---------------------------------------------------------------------------//
/*!
 * \file DTK_ConsistentInterpolation.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation mapping declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATION_HPP
#define DTK_CONSISTENTINTERPOLATION_HPP

#include <DTK_Exception.hpp>
#include <DTK_Rendezvous.hpp>


#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Distributor.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Mesh, class CoordinateField>
ConsistentInterpolation<Mesh,CoordinateField>::ConsistentInterpolation( 
    const RCP_Comm& comm )
    : d_comm( comm )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh, class CoordinateField>
ConsistentInterpolation<Mesh,CoordinateField>::~ConsistentInterpolation()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Setup for interpolation.
 */
template<class Mesh, class CoordinateField>
void ConsistentInterpolation<Mesh,CoordinateField>::setup( 
    const Mesh& mesh, const CoordinateField& coordinate_field )
{
    // Get the global bounding box for the coordinate field.

    // Get the global bounding box for the mesh.

    // Intersect the boxes to get the rendezvous bounding box.

    // Build a rendezvous decomposition.

    // Compute a unique global ordinal for each point in the coordinate field.

    // Build the import map from the global ordinals.

    // Determine the destination of each point in the coordinate field.

    // Via an inverse communication operation, send the global point ordinals
    // to the rendezvous decomposition.

    // Search the rendezvous decomposition with the points.

    // Send the elements / coordinate pairs and the coordinate global ordinals
    // to the local decomposition via an inverse communication operation.

    // Build the data export map from the coordinate ordinals.
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the interpolation.
 */
template<class Mesh, class CoordinateField>
template<class SourceField, class TargetField>
void ConsistentInterpolation<Mesh,CoordinateField>::apply( 
    const FieldEvaluator<Mesh,SourceField>& source_evaluator,
    TargetField& target_space )
{

}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CONSISTENTINTERPOLATION_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolation_def.hpp
//---------------------------------------------------------------------------//

