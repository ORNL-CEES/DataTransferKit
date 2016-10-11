//---------------------------------------------------------------------------//
/*!
 * \file DamperEvaluator.hpp
 * \author Stuart R. Slattery
 * \brief Damper code evaluators for DTK.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EXAMPLE_DAMPEREVALUATOR_HPP
#define DTK_EXAMPLE_DAMPEREVALUATOR_HPP

#include "Damper.hpp"

#include <DTK_FieldContainer.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>

//---------------------------------------------------------------------------//
/*!
  \class DamperEvaluator
  \brief Field Evaluator for the Damper code.
 */
//---------------------------------------------------------------------------//
class DamperEvaluator : public DataTransferKit::FieldEvaluator<
                            int, DataTransferKit::FieldContainer<double>>
{
  public:
    // Typedefs.
    typedef Teuchos::RCP<Damper> RCP_Damper;
    typedef DataTransferKit::MeshContainer<int> mesh_type;
    typedef DataTransferKit::FieldContainer<double> field_type;

    // Constructor.
    DamperEvaluator( const RCP_Damper &damper );

    // Destructor.
    ~DamperEvaluator();

    // Function evaluator.
    field_type evaluate( const Teuchos::ArrayRCP<int> &elements,
                         const Teuchos::ArrayRCP<double> &coords );

  private:
    // Damper instance.
    RCP_Damper d_damper;
};

#endif // end DTK_EXAMPLE_DAMPEREVALUATOR_HPP

//---------------------------------------------------------------------------//
// end DamperEvaluator.hpp
//---------------------------------------------------------------------------//
