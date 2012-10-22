//---------------------------------------------------------------------------//
/*!
 * \file DamperAdapter.hpp
 * \author Stuart R. Slattery
 * \brief Damper code adapters for DTK.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EXAMPLE_DAMPERADAPTER_HPP
#define DTK_EXAMPLE_DAMPERADAPTER_HPP

#include "Damper.hpp"
#include "DamperEvaluator.hpp"

#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTraitsFieldAdapter.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldContainer.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

//---------------------------------------------------------------------------//
/*!
  \class DamperAdapter
  \brief DTK Adpaters for the Damper code.

  This stateless class provides a mechanism for generating mesh and fields
  with the appropriate traits from Damper code data. In addition, this class
  also generates the function evaluator for the damper field.
 */
//---------------------------------------------------------------------------//
class DamperAdapter
{
  public:

    //@{
    //! Typedefs.
    typedef DataTransferKit::MeshContainer<int>                        MeshType;
    typedef DataTransferKit::FieldContainer<double>                    FieldType;
    typedef DataTransferKit::MeshTraits<MeshType>                      MT;
    typedef DataTransferKit::MeshTraits<MeshType>::global_ordinal_type GlobalOrdinal;
    typedef DataTransferKit::FieldEvaluator<GlobalOrdinal,FieldType>   EvaluatorType;
    typedef Teuchos::RCP<EvaluatorType>                                RCP_Evaluator;
    typedef Teuchos::RCP<Damper>                                       RCP_Damper;
    //@}

    // Empty Constructor.
    DamperAdapter()
    { /* ... */ }

    // Destructor.
    ~DamperAdapter()
    { /* ... */ }

    // Get the damper mesh.
    static Teuchos::RCP<DataTransferKit::MeshManager<MeshType> >
    getMesh( const RCP_Damper& damper );

    // Get the damper field evaluator.
    static RCP_Evaluator getFieldEvaluator( const RCP_Damper& damper );

    // Get the damper target coordinates directly from the mesh.
    static Teuchos::RCP<DataTransferKit::FieldManager<MeshType> >
    getTargetCoords( const RCP_Damper& damper );

    // Get the damper target space.
    static Teuchos::RCP<DataTransferKit::FieldManager<FieldType> >
    getTargetSpace( const RCP_Damper& damper );
};

#endif // end DTK_EXAMPLE_DAMPERADAPTER_HPP

//---------------------------------------------------------------------------//
// end DamperAdapter.hpp
//---------------------------------------------------------------------------//

