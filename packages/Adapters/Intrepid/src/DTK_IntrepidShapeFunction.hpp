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
 * \brief DTK_IntrepidShapeFunction.hpp
 * \author Stuart R. Slattery
 * \brief  shape function implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTREPIDSHAPEFUNCTION
#define DTK_INTREPIDSHAPEFUNCTION

#include <unordered_map>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <Intrepid_Basis.hpp>
#include <Intrepid_FieldContainer.hpp>

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class IntrepidShapeFunction
  \brief Intrepid shape function.
*/
//---------------------------------------------------------------------------//
class IntrepidShapeFunction
{
  public:

    /*!
     * \brief Given an topology and a reference point, evaluate the shape
     * function of the topology at that point.
     * \param topology Evaluate the shape function of this topology.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param values Topology shape function evaluated at the reference
     * point. 
     */
    void evaluateValue( 
	const shards::CellTopology& topology,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<double> & values ) const;

    /*!
     * \brief Given an topology and a reference point, evaluate the gradient of
     * the shape function of the topology at that point.
     * \param topology Evaluate the shape function of this topology.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param gradients Topology shape function gradients evaluated at the reference
     * point. Return these ordered with respect to those return by
     * getSupportIds() such that gradients[N][D] gives the gradient value of the
     * Nth support location in the Dth spatial dimension.
     */
    void evaluateGradient( 
	const shards::CellTopology& topology,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<Teuchos::Array<double> >& gradients ) const;

  private:

    // Get the basis of a topology.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
    getIntrepidBasis( const shards::CellTopology& topology ) const;
 
  private:

    // Map of already created shape functions.
    mutable
    std::unordered_map<unsigned,
             Teuchos::RCP<
                 Intrepid::Basis<double,Intrepid::FieldContainer<double>
                                 > > > d_basis;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTREPIDSHAPEFUNCTION

//---------------------------------------------------------------------------//
// end DTK_IntrepidShapeFunction.hpp
//---------------------------------------------------------------------------//
