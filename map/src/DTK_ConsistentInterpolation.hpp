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
#include "DTK_FieldEvaluator"
#include <DTK_MeshTraits.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_Map.hpp>

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
    typedef MT::global_ordinal_type                   global_ordinal_type;
    typedef CoordinateField                           coord_field_type;
    typedef FieldTraits<CoordinateField>              CFT;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos:RCP<CommType>                     RCP_Comm;
    typedef Teuchos::TpetraMap<global_ordinal_type>   TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    //!

    // Constructor.
    ConsistentInterpolation( const RCP_Comm& comm );

    // Destructor.
    ~ConsistentInterpolation();

    // Setup for interpolation.
    void setup( const Mesh& mesh, const CoordinateField& coordinate_field );

    // Apply the interpolation.
    template<class SourceField, class TargetField>
    void apply( const FieldEvaluator<Mesh,SourceField>& source_evaluator,
		TargetField& target_space );

  private:

    // Export field map.
    RCP_TpetraMap d_export_map;

    // Import field map.
    RCP_TpetraMap d_import_map;
};

} // end namespace DataTransferKit

#endif // end DTK_CONSISTENTINTERPOLATION_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolation.hpp
//---------------------------------------------------------------------------//

