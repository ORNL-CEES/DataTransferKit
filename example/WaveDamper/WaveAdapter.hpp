//---------------------------------------------------------------------------//
/*!
 * \file WaveAdapter.hpp
 * \author Stuart R. Slattery
 * \brief Wave code adapters for DTK.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EXAMPLE_WAVEADAPTER_HPP
#define DTK_EXAMPLE_WAVEADAPTER_HPP

#include "Wave.hpp"
#include "WaveEvaluator.hpp"

#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTraitsFieldAdapter.hpp>
#include <DTK_FieldManager.hpp>
#include <DTK_FieldContainer.hpp>
#include <DTK_FieldEvaluator.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

//---------------------------------------------------------------------------//
/*!
  \class WaveAdapter
  \brief DTK Adpaters for the Wave code.

  This stateless class provides a mechanism for generating mesh and fields
  with the appropriate traits from Wave code data. In addition, this class
  also generates the function evaluator for the wave field.
 */
//---------------------------------------------------------------------------//
class WaveAdapter
{
  public:

    //@{
    //! Typedefs.
    typedef DataTransferKit::MeshContainer<int>                  MeshType;
    typedef DataTransferKit::FieldContainer<double>              FieldType;
    typedef DataTransferKit::MeshTraits<MeshType>                MT;
    typedef DataTransferKit::FieldEvaluator<MeshType,FieldType>  EvaluatorType;
    typedef Teuchos::RCP<EvaluatorType>                          RCP_Evaluator;
    typedef Teuchos::RCP<Wave>                                   RCP_Wave;
    //@}

    // Empty Constructor.
    WaveAdapter()
    { /* ... */ }

    // Destructor.
    ~WaveAdapter()
    { /* ... */ }

    // Get the wave mesh.
    static Teuchos::RCP<DataTransferKit::MeshManager<MeshType> >
    getMesh( const RCP_Wave& wave );

    // Get the wave field evaluator.
    static RCP_Evaluator getFieldEvaluator( const RCP_Wave& wave );

    // Get the wave target coordinates directly from the mesh.
    static Teuchos::RCP<DataTransferKit::FieldManager<MeshType> >
    getTargetCoords( const RCP_Wave& wave );

    // Get the wave target space.
    static Teuchos::RCP<DataTransferKit::FieldManager<FieldType> >
    getTargetSpace( const RCP_Wave& wave );
};

#endif // end DTK_EXAMPLE_WAVEADAPTER_HPP

//---------------------------------------------------------------------------//
// end WaveAdapter.hpp
//---------------------------------------------------------------------------//

