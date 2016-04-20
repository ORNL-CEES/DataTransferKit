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
 * \file DTK_CloudDomain.hpp
 * \author Stuart R. Slattery
 * \brief Cloud domain declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLOUDDOMAIN_HPP
#define DTK_CLOUDDOMAIN_HPP

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_SerializationTraits.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class CloudDomain
 * \brief Axis-aligned Cartesian cloud domain container.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class CloudDomain
{

  public:

    // Default constructor.
    CloudDomain();

    // Constructor.
    CloudDomain( const double bounds[2*DIM] );

    // Expand the domain by the given radius.
    void expand( const double radius );

    // Determine if a point is in the domain.
    bool pointInDomain( const Teuchos::ArrayView<const double>& coords ) const;

    // Determine if the given domain intersects this domain.
    bool checkForIntersection( const CloudDomain<DIM>& domain ) const;

    // Get the center of the domain.
    Teuchos::Array<double> center() const;

    // Get the boundaries of the domain.
    Teuchos::ArrayView<const double> bounds() const
    { return Teuchos::ArrayView<const double>(&d_bounds[0],2*DIM); }

 private:

    // Domain boundaries (x_min, x_max, y_min, y_max, z_min, z_max).
    double d_bounds[2*DIM];
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Serialization Traits Specialization.
//---------------------------------------------------------------------------//

namespace Teuchos
{

template<typename Ordinal>
class SerializationTraits<Ordinal, DataTransferKit::CloudDomain<1> >
    : public DirectSerializationTraits<Ordinal, DataTransferKit::CloudDomain<1> >
{ /* ... */ };

template<typename Ordinal>
class SerializationTraits<Ordinal, DataTransferKit::CloudDomain<2> >
    : public DirectSerializationTraits<Ordinal, DataTransferKit::CloudDomain<2> >
{ /* ... */ };

template<typename Ordinal>
class SerializationTraits<Ordinal, DataTransferKit::CloudDomain<3> >
    : public DirectSerializationTraits<Ordinal, DataTransferKit::CloudDomain<3> >
{ /* ... */ };

} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end DTK_CLOUDDOMAIN_HPP

//---------------------------------------------------------------------------//
// end DTK_CloudDomain.hpp
//---------------------------------------------------------------------------//

