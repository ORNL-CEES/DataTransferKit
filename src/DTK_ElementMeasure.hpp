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
 * \file DTK_ElementMeasure.hpp
 * \author Stuart R. Slattery
 * \brief Interface definition for element measure.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ELEMENTMEASURE_HPP
#define DTK_ELEMENTMEASURE_HPP

#include "DTK_FieldTraits.hpp"
#include "DTK_MeshTraits.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
  \class ElementMeasure
  \brief Interface definition for element measures.
 
  Given a list of mesh elements return that element's measure (volume for a
  3D element, area for a 2D element, and length for a 1D element).
*/
//---------------------------------------------------------------------------//
template<class Mesh>
class ElementMeasure
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                mesh_type;
    typedef MeshTraits<Mesh>                    MT;
    typedef typename MT::global_ordinal_type    GlobalOrdinal;
    //@}

    //! Constructor.
    ElementMeasure()
    { /* ... */ }

    //! Destructor.
    virtual ~ElementMeasure()
    { /* ... */ }

    /*!
     * \brief Get the measures of the given elements and return the measures in
     * a FieldTraits container.
     *
     * \param elements A vector of locally valid element global ordinals to
     * get the measures for.
     *
     * \return An array containing the element measures (e.g. volume for 3D
     * elements, area for 2D elements, and length for 1D elements). This
     * returned field is required to be of the same length as the elements
     * input vector and is required to be of a single dimension.
     */
    virtual Teuchos::Array<double>
    measure( const Teuchos::ArrayRCP<GlobalOrdinal>& elements ) = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_ELEMENTMEASURE_HPP

//---------------------------------------------------------------------------//
// end DTK_ElementMeasure.hpp
//---------------------------------------------------------------------------//

