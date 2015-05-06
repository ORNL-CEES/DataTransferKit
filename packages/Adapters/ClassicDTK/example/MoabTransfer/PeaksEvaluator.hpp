//---------------------------------------------------------------------------//
/*!
 * \file PeaksEvaluator.hpp
 * \author Stuart R. Slattery
 * \brief Peaks function evaluator declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PEAKSEVALUATOR_EX_HPP
#define DTK_PEAKSEVALUATOR_EX_HPP

#include "MoabMesh.hpp"
#include "ArrayField.hpp"

#include <DTK_FieldEvaluator.hpp>

class PeaksEvaluator 
    : public DataTransferKit::FieldEvaluator<MoabMesh::global_ordinal_type,
					     ArrayField>
{
  public:

    typedef MoabMesh::global_ordinal_type global_ordinal_type;

    PeaksEvaluator( const MoabMesh& mesh )
	: d_mesh( mesh )
    { /* ... */ }

    ~PeaksEvaluator()
    { /* ... */ }

    ArrayField evaluate( const Teuchos::ArrayRCP<global_ordinal_type>& elements,
			 const Teuchos::ArrayRCP<double>& coords );

  private:

    MoabMesh d_mesh;
};

#endif // end DTK_PEAKSEVALUATOR_EX_HPP

//---------------------------------------------------------------------------//
// end PeaksEvaluator.hpp
//---------------------------------------------------------------------------//

