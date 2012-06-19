//---------------------------------------------------------------------------//
/*!
 * \file DTK_ConsistentInterpolation.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation mapping declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATION_HPP
#define DTK_CONSISTENTINTERPOLATION_HPP

#include "DTK_FieldTraits.hpp"
#include "DTK_FieldEvaluator.hpp"
#include <DTK_MeshTraits.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Import.hpp>

namespace DataTransferKit
{

template<class Mesh, class CoordinateField>
class ConsistentInterpolation
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                      mesh_type;
    typedef MeshTraits<Mesh>                          MT;
    typedef typename MT::global_ordinal_type          global_ordinal_type;
    typedef CoordinateField                           coord_field_type;
    typedef FieldTraits<CoordinateField>              CFT;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;
    typedef Tpetra::Map<global_ordinal_type>          TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    typedef Tpetra::Import<global_ordinal_type>       ImportType;
    typedef Teuchos::RCP<ImportType>                  RCP_Import;
    //!@}

    // Constructor.
    ConsistentInterpolation( const RCP_Comm& comm );

    // Destructor.
    ~ConsistentInterpolation();

    // Setup for interpolation.
    void setup( const Mesh& mesh, const CoordinateField& coordinate_field );

    // Apply the interpolation.
    template<class SourceField, class TargetField>
    void apply( const Teuchos::RCP< 
		    FieldEvaluator<Mesh,SourceField> >& source_evaluator,
		TargetField& target_space );

  private:

    // Build the bounding box for the rendezvous decomposition.
    BoundingBox buildRendezvousBox( const Mesh& mesh, 
				    const BoundingBox& mesh_box,
				    const CoordinateField& coordinate_field, 
				    const BoundingBox& coord_box );

    // Compute globally unique ordinals for the points in the coordinate
    // field.
    Teuchos::Array<global_ordinal_type> computePointOrdinals(
	const CoordinateField& coordinate_field );

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Export field map.
    RCP_TpetraMap d_export_map;

    // Import field map.
    RCP_TpetraMap d_import_map;

    // Field data importer.
    RCP_Import d_data_importer;

    // Local source elements.
    Teuchos::Array<global_ordinal_type> d_source_elements;

    // Local target coords.
    Teuchos::Array<double> d_target_coords;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_ConsistentInterpolation_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_CONSISTENTINTERPOLATION_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolation.hpp
//---------------------------------------------------------------------------//

