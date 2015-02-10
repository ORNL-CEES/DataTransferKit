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
 * \brief DTK_PredicateComposition_impl.hpp
 * \author Stuart R. Slattery
 * \brief Abstract iterator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PREDICATECOMPOSITION_IMPL_HPP
#define DTK_PREDICATECOMPOSITION_IMPL_HPP

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Static Members.
// ---------------------------------------------------------------------------//
// Apply an and operation to two predicates to create a new
// predicate.
template<class ValueType>
std::function<bool(ValueType)>
PredicateComposition::And( const std::function<bool(ValueType)>& func_left,
			   const std::function<bool(ValueType)>& func_right )
{
    return std::bind( std::logical_and<bool>(),
		      std::bind(func_left,std::placeholders::_1),
		      std::bind(func_right,std::placeholders::_1) );
}

//---------------------------------------------------------------------------//
// Apply a or operation to two predicates to create a new predicate.
template<class ValueType>
std::function<bool(ValueType)>
PredicateComposition::Or( const std::function<bool(ValueType)>& func_left,
			  const std::function<bool(ValueType)>& func_right )
{
    return std::bind( std::logical_or<bool>(),
		      std::bind(func_left,std::placeholders::_1),
		      std::bind(func_right,std::placeholders::_1) );
}

//---------------------------------------------------------------------------//
// Apply a not operation to two predicates to create a new
// predicate.
template<class ValueType>
std::function<bool(ValueType)>
PredicateComposition::Not( const std::function<bool(ValueType)>& func )
{
    return std::bind( std::logical_not<bool>(),
		      std::bind(func,std::placeholders::_1) );
}

//---------------------------------------------------------------------------//
// Apply an AndNot operation to create a new predicate.
template<class ValueType>
std::function<bool(ValueType)>
PredicateComposition::AndNot( const std::function<bool(ValueType)>& func_left,
			      const std::function<bool(ValueType)>& func_right )
{
    return std::bind( 
	std::logical_and<bool>(),
	std::bind(func_left,std::placeholders::_1),
	std::bind(Not(func_right),std::placeholders::_1) );
}

//---------------------------------------------------------------------------//
// Apply an OrNot operation to create a new predicate.
template<class ValueType>
std::function<bool(ValueType)>
PredicateComposition::OrNot( const std::function<bool(ValueType)>& func_left,
			     const std::function<bool(ValueType)>& func_right )
{
    return std::bind( 
	std::logical_or<bool>(),
	std::bind(func_left,std::placeholders::_1),
	std::bind(Not(func_right),std::placeholders::_1) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_PREDICATECOMPOSITION_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_PredicateComposition_impl.hpp
//---------------------------------------------------------------------------//
