//---------------------------------------------------------------------------//
/*!
 * \file WaveEvaluator.hpp
 * \author Stuart R. Slattery
 * \brief Wave code evaluators for DTK.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EXAMPLE_WAVEEVALUATOR_HPP
#define DTK_EXAMPLE_WAVEEVALUATOR_HPP

#include "Wave.hpp"

#include <DTK_FieldContainer.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>

//---------------------------------------------------------------------------//
/*!
  \class WaveEvaluator
  \brief Field Evaluator for the Wave code.
 */
//---------------------------------------------------------------------------//
class WaveEvaluator : public DataTransferKit::FieldEvaluator<
                          int, DataTransferKit::FieldContainer<double>>
{
  public:
    // Typedefs.
    typedef Teuchos::RCP<Wave> RCP_Wave;
    typedef DataTransferKit::MeshContainer<int> mesh_type;
    typedef DataTransferKit::FieldContainer<double> field_type;

    // Constructor.
    WaveEvaluator( const RCP_Wave &wave );

    // Destructor.
    ~WaveEvaluator();

    // Function evaluator.
    field_type evaluate( const Teuchos::ArrayRCP<int> &elements,
                         const Teuchos::ArrayRCP<double> &coords );

  private:
    // Wave instance.
    RCP_Wave d_wave;
};

#endif // end DTK_EXAMPLE_WAVEEVALUATOR_HPP

//---------------------------------------------------------------------------//
// end WaveEvaluator.hpp
//---------------------------------------------------------------------------//
