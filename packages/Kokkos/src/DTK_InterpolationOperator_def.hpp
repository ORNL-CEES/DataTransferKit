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
 * \brief DTK_InterpolationOperator_def.hpp
 * \brief Interpolation Operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTERPOLATIONOPERATOR_DEF_HPP
#define DTK_INTERPOLATIONOPERATOR_DEF_HPP

#include "DTK_ConfigDefs.hpp"

#include <Intrepid2_FunctionSpaceTools.hpp>
#include <stdexcept>

namespace DataTransferKit
{
template <typename SC, typename LO, typename GO, typename NO>
InterpolationOperator<SC, LO, GO, NO>::InterpolationOperator(
    Teuchos::RCP<Basis<SC, LO, GO, NO>> basis, DynRankView cell_nodes )
    : _basis( basis )
    , _cell_nodes( cell_nodes )
{
}

// points are the points where we want to interpolate the finite element. The
// points should be defined in the physical space.
template <typename SC, typename LO, typename GO, typename NO>
void InterpolationOperator<SC, LO, GO, NO>::apply( DynRankView values,
                                                   DynRankView coefficients,
                                                   DynRankView phys_points )
{
    unsigned int const n_eval_pts = phys_points.extent( 1 );
    unsigned int const space_dim = phys_points.extent( 2 );
    // Transform the points from the physical space to the reference space
    DynRankView ref_points( "ref_points", phys_points.extent( 0 ), n_eval_pts,
                            space_dim );
    _basis->mapToReferenceFrame( ref_points, phys_points, _cell_nodes );

    // Evaluate the value of the basis functions at the evaluation points
    unsigned int const n_fields = _basis->getCardinality();
    DynRankView ref_basis_values( "ref_basis_values", n_fields, n_eval_pts );
    DynRankView cell_ref_points( "cell_ref_points", n_eval_pts, space_dim );
    for ( unsigned int i = 0; i < n_eval_pts; ++i )
        for ( unsigned int j = 0; j < space_dim; ++j )
            cell_ref_points( i, j ) = ref_points( 0, i, j );

    _basis->getValues( ref_basis_values, cell_ref_points );

    // Transform basis values to physical frame
    DynRankView phys_basis_values( "phys_basis_values", 1, n_fields,
                                   n_eval_pts );
    if ( _basis->getEFunctionSpace() ==
         Intrepid2::EFunctionSpace::FUNCTION_SPACE_HGRAD )
    {
        // This only replicates the value to every cell => (F,P) -> (C,F,P)
        Intrepid2::FunctionSpaceTools<execution_space>::HGRADtransformVALUE(
            phys_basis_values, ref_basis_values );
    }
    else if ( _basis->getEFunctionSpace() ==
              Intrepid2::EFunctionSpace::FUNCTION_SPACE_HDIV )
    {
        // TODO
        // int const space_dim = _cell_topology.getDimension();
        // int const n_cells = 1;
        // DynRankView jacobian( "jacobian", n_cells, n_eval_pts, space_dim,
        //                      space_dim );
        // DynRankView jacobian_det( "jacobian_det", n_cells, n_eval_pts );
        // Intrepid2::CellTools<execution_space>::setJacobian(
        //    jacobian, ref_points, _cell_nodes, _cell_topology );
        // Intrepid2::CellTools<execution_space>::setJacobianDet( jacobian_det,
        //                                                       jacobian );
        //  Intrepid2::FunctionSpaceTools<execution_space>::HDIVtransformVALUE(
        //      phys_basis_values, jacobian, jacobian_det, ref_basis_values );
    }
    else if ( _basis->getEFunctionSpace() ==
              Intrepid2::EFunctionSpace::FUNCTION_SPACE_HCURL )
    {
        // TODO missing inverse jacobian
    }
    else
        throw std::runtime_error( "Not implemented" );

    // TODO Apply field signs for HDIV and HCURL

    // Evaluate finite element at evaluation points
    Intrepid2::FunctionSpaceTools<execution_space>::evaluate(
        values, coefficients, phys_basis_values );
}
}

// Explicit Instantiation Macro.
//---------------------------------------------------------------------------//
#define DTK_INTERPOLATIONOPERATOR_INSTANT( SCALAR, LO, GO, NODE )              \
    template class InterpolationOperator<SCALAR, LO, GO, NODE>;

#endif
