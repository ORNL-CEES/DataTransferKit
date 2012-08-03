//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file DTK_FieldEvaluator.hpp
 * \author Stuart R. Slattery
 * \brief Interface definition for function evaluation kernels.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDEVALUATOR_HPP
#define DTK_FIELDEVALUATOR_HPP

#include "DTK_FieldTraits.hpp"
#include <DTK_MeshTraits.hpp>

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class FieldEvaluator
 * \brief Interface definition for function evaluation kernels.
 */
//---------------------------------------------------------------------------//
template<class Mesh, class Field>
class FieldEvaluator
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                mesh_type;
    typedef MeshTraits<Mesh>                    MT;
    typedef typename MT::global_ordinal_type    GlobalOrdinal;
    typedef Field                               field_type;
    typedef FieldTraits<Field>                  FT;
    typedef typename FT::value_type             data_type;
    //@}

    //! Constructor.
    FieldEvaluator()
    { /* ... */ }

    //! Destructor.
    virtual ~FieldEvaluator()
    { /* ... */ }

    /*!
     * \brief Evaluate the function in the given elements at the given
     * coordinates and return the evaluations in a FieldTraits container.
     *
     * \param elements A vector of valid element handles in which to evaluate
     * the field.
     *
     * \param coords A vector of blocked coordinates 
     * ( x0, x1, x2, ... , xN, y0, y1, y2, ... , yN, z0, z1, z2, ... , zN )
     * at which to evaluate the field. Coordinates ( xN, yN, zN ) should be
     * evaluated in the Nth element in the elements vector.
     *
     * \return A FieldTraits container containing the evaluated function
     * values. This returned field is required to be of the same length as the
     * elements input vector. For those coordinates that can't be evaluated in
     * the given element, return 0 in their position. Field data
     * dimensionality and ordering is specified by field traits.
     */
    virtual Field evaluate( 
	const Teuchos::ArrayRCP<GlobalOrdinal>& elements,
	const Teuchos::ArrayRCP<double>& coords ) = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_FIELDEVALUATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldEvaluator.hpp
//---------------------------------------------------------------------------//

