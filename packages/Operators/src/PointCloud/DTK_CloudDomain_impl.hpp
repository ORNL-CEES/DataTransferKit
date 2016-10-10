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
  * \file DTK_CloudDomain_impl.hpp
  * \author Stuart R. Slattery
  * \brief Cloud domain definition.
  */
//---------------------------------------------------------------------------//

#ifndef DTK_CLOUDDOMAIN_IMPL_HPP
#define DTK_CLOUDDOMAIN_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_CloudDomain.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
template<int DIM>
CloudDomain<DIM>::CloudDomain()
{
    Teuchos::Array<double> zero( 2*DIM, 0.0 );
    std::copy( zero.begin(), zero.end(), d_bounds );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param bounds The bounds of the domain.
 */
template<int DIM>
CloudDomain<DIM>::CloudDomain( const double bounds[2*DIM] )
{
    std::copy( bounds, bounds+2*DIM, d_bounds );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expand the domain by the given radius.
 */
template<>
void CloudDomain<1>::expand( const double radius )
{
    DTK_CHECK( radius >= 0.0 );
    d_bounds[0] -= radius;
    d_bounds[1] += radius;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expand the domain by the given radius.
 */
template<>
void CloudDomain<2>::expand( const double radius )
{
    DTK_CHECK( radius >= 0.0 );
    d_bounds[0] -= radius;
    d_bounds[1] += radius;
    d_bounds[2] -= radius;
    d_bounds[3] += radius;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expand the domain by the given radius.
 */
template<>
void CloudDomain<3>::expand( const double radius )
{
    DTK_CHECK( radius >= 0.0 );
    d_bounds[0] -= radius;
    d_bounds[1] += radius;
    d_bounds[2] -= radius;
    d_bounds[3] += radius;
    d_bounds[4] -= radius;
    d_bounds[5] += radius;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point is in the domain.
 *
 * \param coords Cartesian coordinates to check for point inclusion.
 *
 * \return Return true if the point is in the domain, false if not. A point on
 * the domain boundary will return true.
 */
template<>
bool CloudDomain<1>::pointInDomain(
    const Teuchos::ArrayView<const double>& coords ) const
{
    DTK_REQUIRE( coords.size() == 1 );
    DTK_CHECK( d_bounds[0] <= d_bounds[1] );

    return ( coords[0] >= d_bounds[0] && coords[0] <= d_bounds[1] )
        ? true : false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point is in the domain.
 *
 * \param coords Cartesian coordinates to check for point inclusion.
 *
 * \return Return true if the point is in the domain, false if not. A point on
 * the domain boundary will return true.
 */
template<>
bool CloudDomain<2>::pointInDomain(
    const Teuchos::ArrayView<const double>& coords ) const
{
    DTK_REQUIRE( coords.size() == 2 );
    DTK_CHECK( d_bounds[0] <= d_bounds[1] );
    DTK_CHECK( d_bounds[2] <= d_bounds[3] );

    return ( coords[0] >= d_bounds[0] && coords[0] <= d_bounds[1] &&
             coords[1] >= d_bounds[2] && coords[1] <= d_bounds[3] )
        ? true : false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point is in the domain.
 *
 * \param coords Cartesian coordinates to check for point inclusion.
 *
 * \return Return true if the point is in the domain, false if not. A point on
 * the domain boundary will return true.
 */
template<>
bool CloudDomain<3>::pointInDomain(
    const Teuchos::ArrayView<const double>& coords ) const
{
    DTK_REQUIRE( coords.size() == 3 );
    DTK_CHECK( d_bounds[0] <= d_bounds[1] );
    DTK_CHECK( d_bounds[2] <= d_bounds[3] );
    DTK_CHECK( d_bounds[4] <= d_bounds[5] );

    return ( coords[0] >= d_bounds[0] && coords[0] <= d_bounds[1] &&
             coords[1] >= d_bounds[2] && coords[1] <= d_bounds[3] &&
             coords[2] >= d_bounds[4] && coords[2] <= d_bounds[5] )
        ? true : false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if the given domain intersects this domain.
 *
 * \param domain The domain to check for intersection.
 *
 * \return Return true if the there is an intersection. false if not.
 */
template<>
bool CloudDomain<1>::checkForIntersection(
    const CloudDomain<1>& domain ) const
{
    Teuchos::ArrayView<const double> bounds = domain.bounds();

    return !( ( d_bounds[0] > bounds[1] || d_bounds[1] < bounds[0] ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if the given domain intersects this domain.
 *
 * \param domain The domain to check for intersection.
 *
 * \return Return true if the there is an intersection. false if not.
 */
template<>
bool CloudDomain<2>::checkForIntersection(
    const CloudDomain<2>& domain ) const
{
    Teuchos::ArrayView<const double> bounds = domain.bounds();

    return !( ( d_bounds[0] > bounds[1] || d_bounds[1] < bounds[0] ) ||
              ( d_bounds[2] > bounds[3] || d_bounds[3] < bounds[2] ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if the given domain intersects this domain.
 *
 * \param domain The domain to check for intersection.
 *
 * \return Return true if the there is an intersection. false if not.
 */
template<>
bool CloudDomain<3>::checkForIntersection(
    const CloudDomain<3>& domain ) const
{
    Teuchos::ArrayView<const double> bounds = domain.bounds();

    return !( ( d_bounds[0] > bounds[1] || d_bounds[1] < bounds[0] ) ||
              ( d_bounds[2] > bounds[3] || d_bounds[3] < bounds[2] ) ||
              ( d_bounds[4] > bounds[5] || d_bounds[5] < bounds[4] ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the center of the domain.
 */
template<>
Teuchos::Array<double> CloudDomain<1>::center() const
{
    DTK_CHECK( d_bounds[0] <= d_bounds[1] );
    return Teuchos::Array<double>( 1, (d_bounds[1]+d_bounds[0]) / 2.0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the center of the domain.
 */
template<>
Teuchos::Array<double> CloudDomain<2>::center() const
{
    DTK_CHECK( d_bounds[0] <= d_bounds[1] );
    DTK_CHECK( d_bounds[2] <= d_bounds[3] );

    Teuchos::Array<double> center( 2 );
    center[0] = (d_bounds[1]+d_bounds[0]) / 2.0;
    center[1] = (d_bounds[3]+d_bounds[2]) / 2.0;
    return center;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the center of the domain.
 */
template<>
Teuchos::Array<double> CloudDomain<3>::center() const
{
    DTK_CHECK( d_bounds[0] <= d_bounds[1] );
    DTK_CHECK( d_bounds[2] <= d_bounds[3] );
    DTK_CHECK( d_bounds[4] <= d_bounds[5] );

    Teuchos::Array<double> center( 3 );
    center[0] = (d_bounds[1]+d_bounds[0]) / 2.0;
    center[1] = (d_bounds[3]+d_bounds[2]) / 2.0;
    center[2] = (d_bounds[5]+d_bounds[4]) / 2.0;
    return center;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CLOUDDOMAIN_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_CloudDomain.cpp
//---------------------------------------------------------------------------//

