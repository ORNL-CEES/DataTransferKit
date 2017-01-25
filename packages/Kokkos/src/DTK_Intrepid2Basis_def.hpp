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
 * \brief DTK_Intrepid2Basis_def.hpp
 * \brief Intrepid2Basis_def.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTREPID2BASIS_DEF_HPP
#define DTK_INTREPID2BASIS_DEF_HPP

#include "DTK_ConfigDefs.hpp"

#include <Intrepid2_CellTools.hpp>

namespace DataTransferKit
{

template <typename SC, typename LO, typename GO, typename NO>
Intrepid2Basis<SC, LO, GO, NO>::Intrepid2Basis(
    Teuchos::RCP<Intrepid2::Basis<execution_space>> basis,
    Intrepid2::EFunctionSpace function_space,
    shards::CellTopology const &cell_topology )
    : _basis( basis )
    , _function_space( function_space )
    , _cell_topology( cell_topology )
{
}

template <typename SC, typename LO, typename GO, typename NO>
void Intrepid2Basis<SC, LO, GO, NO>::mapToReferenceFrame(
    DynRankView ref_points, DynRankView phys_points, DynRankView cell_nodes )
{
    Intrepid2::CellTools<execution_space>::mapToReferenceFrame(
        ref_points, phys_points, cell_nodes, _cell_topology );
}

template <typename SC, typename LO, typename GO, typename NO>
void Intrepid2Basis<SC, LO, GO, NO>::getValues(
    DynRankView ref_basis_values, DynRankView const cell_ref_points )
{
    _basis->getValues( ref_basis_values, cell_ref_points,
                       Intrepid2::OPERATOR_VALUE );
}

template <typename SC, typename LO, typename GO, typename NO>
unsigned int Intrepid2Basis<SC, LO, GO, NO>::getCardinality()
{
    return _basis->getCardinality();
}

template <typename SC, typename LO, typename GO, typename NO>
Intrepid2::EFunctionSpace Intrepid2Basis<SC, LO, GO, NO>::getEFunctionSpace()
{
    return _function_space;
}
}

// Explicit Instantiation Macro.
//---------------------------------------------------------------------------//
#define DTK_INTREPID2BASIS_INSTANT( SCALAR, LO, GO, NODE )                     \
    template class Intrepid2Basis<SCALAR, LO, GO, NODE>;

#endif
