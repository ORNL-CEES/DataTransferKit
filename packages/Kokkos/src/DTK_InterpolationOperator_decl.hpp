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
 * \brief DTK_InterpolationOperator_decl.hpp
 * \brief Interpolation Operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTERPOLATIONOPERATOR_DECL_HPP
#define DTK_INTERPOLATIONOPERATOR_DECL_HPP

#include "DTK_ConfigDefs.hpp"

#include <Intrepid2_Basis.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>
#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{
template <typename SC, typename LO, typename GO, typename NO>
class InterpolationOperator
{
  public:
    using scalar_type = SC;
    using local_ordinal_type = LO;
    using global_ordinal_type = GO;
    using node_type = NO;
    using device_type = typename NO::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;
    typedef Kokkos::Experimental::DynRankView<double, execution_space>
        DynRankView;

    InterpolationOperator(
        Teuchos::RCP<Intrepid2::Basis<execution_space>> basis,
        shards::CellTopology const &_cell_topology, DynRankView cell_nodes,
        Intrepid2::EFunctionSpace function_space );

    void apply( DynRankView value, DynRankView coefficients,
                DynRankView phys_points );

  private:
    Teuchos::RCP<Intrepid2::Basis<execution_space>> _basis;
    shards::CellTopology _cell_topology;
    DynRankView _cell_nodes;
    Intrepid2::EFunctionSpace _function_space;
};
}

#endif
