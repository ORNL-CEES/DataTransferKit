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
 * \file DTK_EvaluationPoint_impl.hpp
 * \author Stuart R. Slattery
 * \brief EvaluationPoint class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EVALUATIONPOINT_IMPL_HPP
#define DTK_EVALUATIONPOINT_IMPL_HPP

#include <algorithm>

#include "DTK_DBC.hpp"
#include "DTK_Serializer.hpp"

#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Deserializer constructor.
 */
template<class Ordinal>
EvaluationPoint<Ordinal>::EvaluationPoint( const Teuchos::ArrayView<char>& buffer )
{
    DTK_REQUIRE( Teuchos::as<std::size_t>(buffer.size()) == d_packed_bytes );

    Deserializer ds;
    ds.setBuffer( d_packed_bytes, buffer.getRawPtr() );
    ds >> d_gid >> d_coords;

    DTK_ENSURE( ds.getPtr() == ds.end() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack the history into a buffer.
 */
template<class Ordinal>
Teuchos::Array<char> EvaluationPoint<Ordinal>::pack() const
{
    DTK_REQUIRE( d_packed_bytes );
    DTK_REQUIRE( d_packed_bytes > 0 );

    Teuchos::Array<char> buffer( d_packed_bytes );

    Serializer s;
    s.setBuffer( d_packed_bytes, buffer.getRawPtr() );
    s << d_gid << d_coords;

    DTK_ENSURE( s.getPtr() == s.end() );
    return buffer;
}

//---------------------------------------------------------------------------//
// Static members.
//---------------------------------------------------------------------------//
template<class Ordinal>
std::size_t EvaluationPoint<Ordinal>::d_packed_bytes = 0;

//---------------------------------------------------------------------------//
/*!
 * \brief Set the byte size of the packed history state.
 */
template<class Ordinal>
void EvaluationPoint<Ordinal>::setByteSize( std::size_t size_rng_state )
{
    d_packed_bytes = sizeof(Ordinal) + DIM*sizeof(double);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of bytes in the packed history state.
 */
template<class Ordinal>
std::size_t EvaluationPoint<Ordinal>::getPackedBytes()
{
    DTK_REQUIRE( d_packed_bytes );
    return d_packed_bytes;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_EVALUATIONPOINT_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_EvaluationPoint_impl.hpp
//---------------------------------------------------------------------------//

