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
 * \file   DTK_IntrepidSideCell.hpp
 * \author Stuart Slattery
 * \brief  Manager for Intrepid cell-level operations on cell sides.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTREPIDSIDECELL_HPP
#define DTK_INTREPIDSIDECELL_HPP

#include "DTK_IntrepidCell.hpp"

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class IntrepidSideCell
 * \brief Manager for Intrepid cell-level operations on cell sides.
 */
//---------------------------------------------------------------------------//
class IntrepidSideCell : public IntrepidCell
{
  public:
    //@{
    //! Typedefs.
    typedef IntrepidCell Base;
    typedef Base::Scalar Scalar;
    typedef Base::MDArray MDArray;
    //@}

  public:
    // Constructor.
    IntrepidSideCell( const shards::CellTopology &side_topology,
                      const unsigned side_id,
                      const shards::CellTopology &parent_topology,
                      const unsigned degree );

    // Update the cell state of the object for the current cell node
    // coordinates.
    void updateCellState();

    // Given a set of coordinates in the reference frame of the side of the
    // parent cell, map them to the physical frame.
    void mapToCellPhysicalFrame( const MDArray &parametric_coords,
                                 MDArray &physical_coords );

    // Compute the physical normals of the side at the integration points
    // (IP,DIM).
    void getPhysicalSideNormalsAtIntegrationPoints( MDArray &side_normals );

    // Compute the physical normals of the side at a given reference point.
    void
    getPhysicalSideNormalsAtReferencePoint( const MDArray &parametric_coords,
                                            MDArray &side_normals );

  private:
    // Reference cell side id on the parent cell.
    unsigned d_side_id;

    // Parent cell topology.
    shards::CellTopology d_parent_topology;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTREPIDSIDECELL_HPP

//---------------------------------------------------------------------------//
// end DTK_IntrepidSideCell.hpp
//---------------------------------------------------------------------------//
