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
 * \file DTK_Classic_FieldContainer.hpp
 * \author Stuart R. Slattery
 * \brief A simple default field container with a field traits
 * implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_Classic_FIELDCONTAINER_HPP
#define DTK_Classic_FIELDCONTAINER_HPP

#include "DTK_Classic_FieldTraits.hpp"

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{
namespace Classic
{
//---------------------------------------------------------------------------//
/*!
 * \class FieldContainer
 * \brief A default field implementation.

 The FieldContainer is a default implementation of FieldTraits for clients. Its
 data mimics the structure of FieldTraits with constructor arguments expected
 to have the data layout of a FieldTraits object. 
 */
//---------------------------------------------------------------------------//
template<typename Scalar>
class FieldContainer
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::ArrayRCP<Scalar>                 field_type;
    typedef typename field_type::value_type           value_type;
    typedef typename field_type::size_type            size_type;
    typedef typename field_type::iterator             iterator;
    typedef typename field_type::const_iterator       const_iterator;

    //@}
    
    //! Default Constructor.
    FieldContainer()
    { /* ... */ }

    /*!
     * \brief Constructor.
     *
     * \param data The data for the field, blocked by dimension.
     *
     * \param dimension The dimension of the field.
     */
    FieldContainer( const Teuchos::ArrayRCP<Scalar>& data,
		    const int dimension )
	: d_data( data )
	, d_dimension( dimension )
    { /* ... */ }

    //! Get the dimension of the field.
    int dim() const
    { return d_dimension; }

    //! Get the size of the field.
    size_type size() const
    { return d_data.size(); }

    //! Returns if the field is empty.
    bool empty() const
    { 
	if ( 0 == d_data.size() ) return true;
	else return false; 
    }

    //@{
    //! Get the beginning of the field.
    iterator begin()
    { return d_data.begin(); }

    const_iterator begin() const
    { return d_data.begin(); }
    //@}

    //@{
    //! Get the end of the field.
    iterator end()
    { return d_data.end(); }

    const_iterator end() const
    { return d_data.end(); }
    //@}

    //! Get the data.
    Teuchos::ArrayRCP<Scalar> getData()
    { return d_data; }

 private:

    // Field data.
    Teuchos::ArrayRCP<Scalar> d_data;

    // Field dimension;
    int d_dimension;
};

//---------------------------------------------------------------------------//
/*
 * \brief FieldTraits specialization for the field container.
 */
//---------------------------------------------------------------------------//
template<typename Scalar>
class FieldTraits< FieldContainer<Scalar> >
{
  public:

    typedef FieldContainer<Scalar>                         Container;
    typedef typename Container::field_type                 field_type;
    typedef typename Container::value_type                 value_type;
    typedef typename Container::size_type                  size_type;
    typedef typename Container::iterator                   iterator;
    typedef typename Container::const_iterator             const_iterator;

    static inline int dim( const Container& container )
    { return container.dim(); }

    static inline size_type size( const Container& container )
    { return container.size(); }

    static inline bool empty( const Container& container )
    { return container.empty(); }

    static inline iterator begin( Container& container )
    { return container.begin(); }

    static inline const_iterator begin( const Container& container )
    { return container.begin(); }

    static inline iterator end( Container& container )
    { return container.end(); }

    static inline const_iterator end( const Container& container )
    { return container.end(); }
};

//---------------------------------------------------------------------------//

} // end namespace Classic
} // end namespace DataTransferKit

#endif // end DTK_Classic_FIELDCONTAINER_HPP

//---------------------------------------------------------------------------//
// end DTK_Classic_FieldContainer.hpp
//---------------------------------------------------------------------------//
