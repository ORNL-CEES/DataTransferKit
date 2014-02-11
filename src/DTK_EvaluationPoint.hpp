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
 * \file DTK_EvaluationPoint.hpp
 * \author Stuart R. Slattery
 * \brief Evaluation point class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EVALUATIONPOINT_HPP
#define DTK_EVALUATIONPOINT_HPP

#include "DTK_BufferDataTraits.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Evaluation point
 * \brief Encapsulation of an evaluation point in a parallel calculation
 */
//---------------------------------------------------------------------------//
template<class Ordinal, int DIM>
class EvaluationPoint
{
  public:

    //@{
    //! Typedefs.
    typedef Ordinal                                   ordinal_type;
    typedef RNGControl::RNG                           RNG;
    //@}

    //! Default constructor.
    EvaluationPoint()
	: d_gid( Teuchos::OrdinalTraits<Ordinal>::invalid() )
    { /* ... */ }

    //! State constructor.
    EvaluationPoint( const Ordinal gid, const double coords[DIM] )
	: d_gid( gid )
    {
	std::copy( &coords[0], &coords[DIM]+1, &d_coords[0] );
    }

    // Deserializer constructor.
    explicit EvaluationPoint( const Teuchos::ArrayView<char>& buffer );

    // Destructor.
    ~EvaluationPoint()
    { /* ... */ }

    // Pack the evaluation point into a buffer.
    Teuchos::Array<char> pack() const;

    //! Set the evaluation point gid.
    inline void setGid( const Ordinal gid )
    { d_gid = gid; }

    //! Get the evaluation point gid.
    inline Ordinal gid() const 
    { return d_gid; }

    //! Set the evaluation point coordinates.
    inline void setCoords( const double coords[DIM] )
    { d_coords = coords; }

    //! Get the evaluation point coords.
    inline Teuchos::ArrayView<double> coords() const 
    { return Teuchos::ArrayView<double>(d_coords,DIM); }

    //! Get the dimension of the point.
    inline int dim() const
    { return DIM; }

  public:

    // Set the byte size of the packed evaluation_point state.
    static void setByteSize();

    // Get the number of bytes in the packed evaluation_point state.
    static std::size_t getPackedBytes();

  private:

    // Point global id.
    Ordinal d_gid;

    // Point coordinates.
    double d_coords[DIM];

  private:

    // Packed size of point in bytes.
    static std::size_t d_packed_bytes;
};

//---------------------------------------------------------------------------//
// BufferDataTraits Implementation.
//---------------------------------------------------------------------------//
template<class Ordinal>
class BufferDataTraits<EvaluationPoint<Ordinal> >
{
  public:

    //@{
    //! Typedefs.
    typedef EvaluationPoint<Ordinal>                     data_type;
    typedef typename data_type::ordinal_type             ordinal_type;
    //@}

    /*!
     * \brief Create a point from a buffer.
     */
    static Teuchos::RCP<data_type> 
    createFromBuffer( const Teuchos::ArrayView<char>& buffer )
    { 
	return Teuchos::rcp( new data_type(buffer) );
    }

    /*!
     * \brief Pack the point into a buffer.
     */
    static Teuchos::Array<char> pack( const data_type& point )
    {
	return point.pack();
    }

    /*!
     * \brief Set the byte size of the packed point state.
     */
    static void setByteSize( std::size_t size_rng_state )
    {
	data_type::setByteSize( size_rng_state );
    }

    /*!
     * \brief Get the number of bytes in the packed point state.
     */
    static std::size_t getPackedBytes()
    {
	return data_type::getPackedBytes();
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_EvaluationPoint_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_EVALUATIONPOINT_HPP

//---------------------------------------------------------------------------//
// end DTK_EvaluationPoint.hpp
//---------------------------------------------------------------------------//

