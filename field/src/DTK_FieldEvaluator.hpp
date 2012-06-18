//---------------------------------------------------------------------------//
/*!
 * \file DTK_FieldEvaluator.hpp
 * \author Stuart R. Slattery
 * \brief Protocol definition for field evaluation kernels.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDEVALUATOR_HPP
#define DTK_FIELDEVALUATOR_HPP

#include "DTK_FieldTraits.hpp"
#include <DTK_MeshTraits.hpp>

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

template<class Mesh, class DataField>
class FieldEvaluator
{
  public:

    //@{
    //! Typedefs.
    typedef MeshTraits<Mesh>                    MT;
    typedef typename MT::global_ordinal_type    global_ordinal_type;
    typedef DataField                           data_field_type;
    typedef FieldTraits<DataField>              FT;
    typedef typename FT::value_type             data_type;
    //@}

    //! Constructor.
    FieldEvaluator()
    { /* ... */ }

    //! Destructor.
    virtual ~FieldEvaluator()
    { /* ... */ }

    /*!
     * \brief Evaluate the field in the given elements at the given
     * coordinates and return the evaluations in a DataField.
     * \param elements A vector of element handles in which to evaluate
     * the field.
     * \param coords A vector of blocked coordinates 
     * ( x0, x1, x2, ... , xN, y0, y1, y2, ... , yN, z0, z1, z2, ... , zN )
     * at which to evaluate the field. Coordinates ( xN, yN, zN ) should be
     * evaluated in the Nth element in the elements vector.
     * \return Return a DataField containing the evaluated field
     * values. This returned field is required to be of the same length as
     * the elements input vector. For those coordinates that can't be
     * evaluated in the given element, return 0 in their position. Field data
     * dimensionality and ordering is specified by field traits.
     */
    virtual DataField evaluate( 
	const Teuchos::ArrayRCP<global_ordinal_type>& elements,
	const Teuchos::ArrayRCP<double>& coords ) = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_FIELDEVALUATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldEvaluator.hpp
//---------------------------------------------------------------------------//

