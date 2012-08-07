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
 * \file DTK_FieldTraits.hpp
 * \author Stuart R. Slattery
 * \brief Traits declaration for field types.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDTRAITS_HPP
#define DTK_FIELDTRAITS_HPP

#include <iterator>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
//---------------------------------------------------------------------------//
template<typename UndefinedFieldType>
struct UndefinedFieldTraits
{
    static inline UndefinedFieldType notDefined() 
    { return UndefinedFieldType::this_type_is_missing_a_specialization(); }
};

//---------------------------------------------------------------------------//
/*!
 * \class FieldTraits
 * \brief Field traits definitions.
  
 These traits correlate to the basic concept of a field within DTK. A field
 can contain anything in an array, but it must store its objects in contiguous
 memory that is blocked by dimension and accessible by random access
 iterators.

 In the most general sense, a field refers to the degrees of freedom computed
 by a physics code or the responses derived from those degrees of freedom that
 have been discretized across the domain. The field is implicitly bound to the
 geometric domain through the degrees of freedom and their association with a
 mesh or other geometric components. In a physics simulation, examples of
 degrees of freedom include pressure and velocity distributions and examples
 of computed reponses include heat flux or reaction rates. In order to access
 DTK field services, a subset of information needed to describe the field is
 required. A field has a dimension of arbitrary size. As examples, for scalar
 fields this dimension is 1, for 3-vectors (such as the velocity example above
 in a 3 dimensional computation) the dimension is 3, and for a 3x3 tensor the
 dimension is 9.  All local instances of the field must have the same
 dimension. A field can have an arbitrary number of local degrees of freedom
 and this size can differ from local domain to local domain. No knowledge of
 the global field decomposition is required, however, it must exist on a
 single communicator.
 */
//---------------------------------------------------------------------------//
template<typename FieldType>
class FieldTraits
{
  public:

    //! Typedef for field type.
    typedef FieldType field_type;

    //! Typedef for value type. The field type must implement
    //! Teuchos::ScalarTraits.
    typedef typename FieldType::value_type value_type;

    //! Typedef for size type.
    typedef typename FieldType::size_type size_type;
    
    //{@
    //! Typedef for field value random access iterators.
    typedef typename 
    std::iterator<std::random_access_iterator_tag,value_type> 
    iterator;

    typedef typename 
    std::iterator<std::random_access_iterator_tag,const value_type> 
    const_iterator;
    //@}

    /*! 
     * \brief Returns the dimensionality of the field ( i.e. 1 for a scalar, 3
     * for a 3 vector, 9 for a 3x3 tensor, etc. ).
     */
    static inline int dim( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    /*! 
     * \brief Returns the total number of elements in the field.
     */
    static inline size_type size( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    /*! 
     * \brief Returns if the field is empty.
     */
    static inline bool empty( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    //@{
    /*! 
     * \brief Returns the iterator to the beginning of the field. The data is
     * required to be blocked by dimensions { d0, d1, d2, ... , dM } as
     * { d0_0, d0_1, ... , d0_N, d1_0, d1_1, ... , d1_N, ... , dM_N }
     */
    static inline iterator begin( FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    static inline const_iterator begin( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }
    //@}

    //@{
    /*! 
     * \brief Returns the iterator to the end of the field. The data is
     * required to be blocked by dimensions { d0, d1, d2, ... , dM } as
     * { d0_0, d0_1, ... , d0_N, d1_0, d1_1, ... , d1_N, ... , dM_N }
     */
    static inline iterator end( FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    static inline const_iterator end( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }
    //@}
};

} // end namespace DataTransferKit

#endif // end DTK_FIELDTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldTraits.hpp
//---------------------------------------------------------------------------//

