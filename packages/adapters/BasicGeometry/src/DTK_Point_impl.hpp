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
 * \file DTK_Point_impl.hpp
 * \author Stuart R. Slattery
 * \brief Point definition
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POINT_IMPL_HPP
#define DTK_POINT_IMPL_HPP

#include <limits>

#include "DTK_DBC.hpp"
#include "DTK_DataSerializer.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
template<int DIM>
Point<DIM>::Point()
{
    d_point_impl = Teuchos::rcp( new PointImpl<DIM>() );
    this->b_entity_impl = d_point_impl;
}

//---------------------------------------------------------------------------//
// Array constructor.
template<int DIM>
Point<DIM>::Point( const EntityId global_id, 
	      const int owner_rank,
	      const Teuchos::Array<double>& coordinates )
{
    d_point_impl = 
	Teuchos::rcp( new PointImpl<DIM>(global_id,owner_rank,coordinates) );
    this->b_entity_impl = d_point_impl;
}

//---------------------------------------------------------------------------//
// Destructor.
template<int DIM>
Point<DIM>::~Point()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the coordinates of the point.
template<int DIM>
void Point<DIM>::getCoordinates( 
    Teuchos::ArrayView<const double>& coordinates ) const
{ 
    d_point_impl->getCoordinates( coordinates );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the point description to an ostream.
 *
 * \return The ostream.
 */
template<int DIM>
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Point<DIM>& p)
{
    Teuchos::ArrayView<const double> coords;
    p.getCoordinates( coords );
    os << "Point: d_global_id=" << p.id()
       << ",d_owner_rank=" << p.ownerRank()
       << ",d_coordinates=" << coords;

  return os;
}

//---------------------------------------------------------------------------//
// Static Members.
//---------------------------------------------------------------------------//
// Get the byte size of the point.
template<int DIM>
std::size_t Point<DIM>::byteSize()
{
    return PointImpl<DIM>::byteSize();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_POINT_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_Point_impl.hpp
//---------------------------------------------------------------------------//

