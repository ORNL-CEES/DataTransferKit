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
 * \file   DTK_IntrepidBasisFactory.cpp
 * \author Stuart Slattery
 * \brief  Factory for Intrepid basis functions.
 */
//---------------------------------------------------------------------------//

#include "DTK_IntrepidBasisFactory.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_RCP.hpp>

#include <Shards_BasicTopologies.hpp>

#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid_HGRAD_PYR_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HGRAD_WEDGE_C1_FEM.hpp"
#include "Intrepid_HGRAD_WEDGE_C2_FEM.hpp"
#include <Intrepid_Basis.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class IntrepidBasisFactory
 * \brief Factory for Intrepid basis functions.
 */
//---------------------------------------------------------------------------//
Teuchos::RCP<Intrepid::Basis<double, Intrepid::FieldContainer<double>>>
IntrepidBasisFactory::create( const shards::CellTopology &cell_topo )
{

    Teuchos::RCP<Intrepid::Basis<double, Intrepid::FieldContainer<double>>>
        basis;

    switch ( cell_topo.getKey() )
    {

    case shards::Line<2>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_LINE_C1_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Triangle<3>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TRI_C1_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Triangle<6>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TRI_C2_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Quadrilateral<4>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Quadrilateral<9>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Tetrahedron<4>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TET_C1_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Tetrahedron<10>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TET_C2_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Hexahedron<8>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_HEX_C1_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Hexahedron<27>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_HEX_C2_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Wedge<6>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_WEDGE_C1_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Wedge<18>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_WEDGE_C2_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    case shards::Pyramid<5>::key:
    case shards::Pyramid<13>::key:
    case shards::Pyramid<14>::key:
        basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_PYR_C1_FEM<
                              double, Intrepid::FieldContainer<double>>() );
        break;

    default:
        bool topology_supported = false;
        DTK_INSIST( topology_supported );
        break;
    }

    return basis;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_IntrepidBasisFactory.cpp
//---------------------------------------------------------------------------//
