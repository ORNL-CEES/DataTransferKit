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

namespace DataTransferKit
{

template<class Mesh, class DataField>
class FieldEvaluator
{
  public:

    //@{
    //! Typedefs.
    typedef MeshTraits<Mesh>                 MT;
    typedef typename MT::handle_type         handle_type;
    typedef DataField                        data_field_type;
    typedef FieldTraits<DataField>           FT;
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
     * \param coords A vector of interleaved coordinates 
     * ( x0, y0, z0, ..., xN, yN, zN ) at which to evaluate the
     * field. Coordinates ( x_n, y_n, z_n ) should be evaluated in the nth
     * element in the elements vector.
     * \return Return a DataField containing the evaluated field
     * values. This returned field is expected to be of the same length as
     * the elements input vector. For those coordinates that can't be
     * evaluated in the given element, return 0 in their position.
     */
    virtual DataField evaluate( const std::vector<handle_type>& elements,
				const std::vector<double>& coords ) = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_FIELDEVALUATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldEvaluator.hpp
//---------------------------------------------------------------------------//

