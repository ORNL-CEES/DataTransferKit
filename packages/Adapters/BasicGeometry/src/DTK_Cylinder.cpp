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
  * \file DTK_Cylinder.cpp
  * \author Stuart R. Slattery
  * \brief Bounding box definition.
  */
//---------------------------------------------------------------------------//

#include <cmath>

#include "DTK_Cylinder.hpp"
#include "DTK_CylinderImpl.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
Cylinder::Cylinder()
{
    this->b_entity_impl = Teuchos::rcp( new CylinderImpl() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param length Length of cylinder along Z-axis.
 *
 * \param radius Radius of cylinder.
 *
 * \param centroid_x Centroid X-coordinate.
 *
 * \param centroid_y Centroid Y-coordinate.
 *
 * \param centroid_z Centroid Z-coordinate.
 */
Cylinder::Cylinder( const EntityId global_id, 
		    const int owner_rank, 
		    const int block_id,
		    const double length, 
		    const double radius,
		    const double centroid_x, 
		    const double centroid_y, 
		    const double centroid_z )
{
    this->b_entity_impl = Teuchos::rcp( 
	new CylinderImpl(global_id,owner_rank,block_id,length,radius,
			 centroid_x,centroid_y,centroid_z) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Cylinder::~Cylinder()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the length of the cylinder.
double Cylinder::length() const
{ 
    return 
	Teuchos::rcp_dynamic_cast<CylinderImpl>(this->b_entity_impl)->length();
}

//---------------------------------------------------------------------------//
// Get the radius of the cylinder.
double Cylinder::radius() const
{ 
    return 
	Teuchos::rcp_dynamic_cast<CylinderImpl>(this->b_entity_impl)->radius();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the measure of the box.
 *
 * \return Return the measure of the box.
 */
double Cylinder::measure() const
{
    return Teuchos::rcp_dynamic_cast<CylinderImpl>(this->b_entity_impl)->measure();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the box.
 *
 * \return The centroid coordinates.
 */
void Cylinder::centroid( const Teuchos::ArrayView<double>& centroid ) const
{
    Teuchos::rcp_dynamic_cast<CylinderImpl>(this->b_entity_impl)->centroid(centroid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a point to the reference space of an entity. Return the
 */
bool Cylinder::mapToReferenceFrame( 
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    return Teuchos::rcp_dynamic_cast<CylinderImpl>(
	this->b_entity_impl)->mapToReferenceFrame(point,reference_point);
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Determine if a reference point is in the parameterized space of
 * an entity.
 */
bool Cylinder::checkPointInclusion( 
    const double tolerance,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    return Teuchos::rcp_dynamic_cast<CylinderImpl>(
	this->b_entity_impl)->checkPointInclusion(tolerance,reference_point);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a reference point to the physical space of an entity.
 */
void Cylinder::mapToPhysicalFrame( 
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    Teuchos::rcp_dynamic_cast<CylinderImpl>(
	this->b_entity_impl)->mapToPhysicalFrame(reference_point,point);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the cylinder description to an ostream.
 *
 * \return The ostream.
 */
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Cylinder& c)
{
    Teuchos::Array<double> centroid(3);
    c.centroid( centroid() );
    os << "Cylinder: length=" << c.length() << ",radius=" << c.radius()
       << ", centroid=(" << centroid[0] << "," << centroid[1]
       << "," << centroid[2] << ")";

  return os;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Cylinder.cpp
//---------------------------------------------------------------------------//

