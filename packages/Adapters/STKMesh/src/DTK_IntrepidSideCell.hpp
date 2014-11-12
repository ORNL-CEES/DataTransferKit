//---------------------------------------------------------------------------//
/*
  Copyright (c) 2013, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
template<typename MDArray>
class IntrepidSideCell : public IntrepidCell<MDArray>
{
  public:

    //@{
    //! Typedefs.
    typedef IntrepidCell<MDArray> Base;
    typedef typename Base::Scalar Scalar;
    //@}

  public:

    // Constructor.
    IntrepidSideCell( const shards::CellTopology& side_topology,
		      const unsigned side_id,
		      const shards::CellTopology& parent_topology,
		      const unsigned degree );

    // Destructor.
    ~IntrepidSideCell();

    // Update the cell state of the object for the current cell node
    // coordinates.
    void updateCellState();

    // Given a set of coordinates in the reference frame of the side of the
    // parent cell, map them to the physical frame.
    void mapToCellPhysicalFrame( const MDArray& parametric_coords,
				 MDArray& physical_coords );

    // Compute the physical normals of the side.
    void getPhysicalSideNormals( MDArray& side_normals );

  private:

    // Reference cell side id over which cell is occuring.
    unsigned d_side_id;

    // Parent cell topology.
    shards::CellTopology d_parent_topology;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_IntrepidSideCell_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_INTREPIDSIDECELL_HPP

//---------------------------------------------------------------------------//
// end DTK_IntrepidSideCell.hpp
//---------------------------------------------------------------------------//

