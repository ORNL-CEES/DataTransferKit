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
 * \brief DTK_PredicateComposition.hpp
 * \author Stuart R. Slattery
 * \brief Tools for predicate composition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PREDICATECOMPOSITION_HPP
#define DTK_PREDICATECOMPOSITION_HPP

#include <functional>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class PredicateComposition
  \brief Tools for predicate composition.

  A stateless class of tools for predicate composition.
*/
//---------------------------------------------------------------------------//
class PredicateComposition
{
  public:

    //! Predicate alias.
    template<class T>
    using Predicate = std::function<bool(T)>;

    /*!
     * \brief Constructor.
     */
    PredicateComposition() { /* ... */ }

    // Apply an and operation to two predicates to create a new
    // predicate.
    template<class ValueType>
    static Predicate<ValueType>
    And( const Predicate<ValueType>& func_left,
         const Predicate<ValueType>& func_right );

    // Apply an or operation to two predicates to create a new predicate.
    template<class ValueType>
    static Predicate<ValueType>
    Or( const Predicate<ValueType>& func_left,
        const Predicate<ValueType>& func_right );

    // Apply a not operation to a predicate to create a new
    // predicate.
    template<class ValueType>
    static Predicate<ValueType>
    Not( const Predicate<ValueType>& func );

    // Apply an AndNot operation to create a new predicate.
    template<class ValueType>
    static Predicate<ValueType>
    AndNot( const Predicate<ValueType>& func_left,
            const Predicate<ValueType>& func_right );

    // Apply an OrNot operation to create a new predicate.
    template<class ValueType>
    static Predicate<ValueType>
    OrNot( const Predicate<ValueType>& func_left,
           const Predicate<ValueType>& func_right );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_PredicateComposition_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_PREDICATECOMPOSITION_HPP

//---------------------------------------------------------------------------//
// end DTK_PredicateComposition.hpp
//---------------------------------------------------------------------------//
