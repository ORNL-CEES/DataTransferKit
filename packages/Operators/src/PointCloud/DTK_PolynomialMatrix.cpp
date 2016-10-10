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
 * \file   DTK_PolynomialMatrix.cpp
 * \author Stuart R. Slattery
 * \brief  Polynomial matrix.
 */
//---------------------------------------------------------------------------//

#include "DTK_DBC.hpp"
#include "DTK_PolynomialMatrix.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <mpi.h>
#endif

#include <Tpetra_Export.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
PolynomialMatrix::PolynomialMatrix(
    const Teuchos::RCP<const Tpetra::MultiVector<double,int,SupportId> >& polynomial,
    const Teuchos::RCP<const Tpetra::Map<int,SupportId> >& domain_map,
    const Teuchos::RCP<const Tpetra::Map<int,SupportId> >& range_map )
    : d_comm( polynomial->getMap()->getComm() )
    , d_polynomial( polynomial )
    , d_domain_map( domain_map )
    , d_range_map( range_map )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Apply operation.
void PolynomialMatrix::apply(
    const Tpetra::MultiVector<double,int,SupportId> &X,
    Tpetra::MultiVector<double,int,SupportId> &Y,
    Teuchos::ETransp mode,
    double alpha,
    double beta ) const
{
    DTK_REQUIRE( d_domain_map->isSameAs(*(X.getMap())) );
    DTK_REQUIRE( d_range_map->isSameAs(*(Y.getMap())) );
    DTK_REQUIRE( X.getNumVectors() == Y.getNumVectors() );

    // Get the size of the problem and view of the local vectors.
    int poly_size = d_polynomial->getNumVectors();
    int num_vec = X.getNumVectors();
    int local_length = d_polynomial->getLocalLength();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > poly_view =
        d_polynomial->get2dView();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > y_view =
        Y.get2dViewNonConst();

    // Scale Y by beta.
    Y.scale( beta );

    // No transpose.
    if ( Teuchos::NO_TRANS == mode )
    {
        // Broadcast the polynomial components of X from the root rank.
        Teuchos::Array<double> x_poly( poly_size * num_vec, 0.0 );
        if ( 0 == d_comm()->getRank() )
        {
            Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > x_view =
                X.get2dView();
            for ( int n = 0; n < num_vec; ++n )
            {
                std::copy( &x_view[n][0], &x_view[n][0] + poly_size,
                           &x_poly[n*poly_size] );
            }
        }
        Teuchos::broadcast( *d_comm, 0, x_poly() );

        // Do the local mat-vec.
        int stride = 0;
        for ( int n = 0; n < num_vec; ++n )
        {
            stride = n*poly_size;

            for ( int i = 0; i < local_length; ++i )
            {
                for ( int p = 0; p < poly_size; ++p )
                {
                    y_view[n][i] += alpha * poly_view[p][i] * x_poly[stride + p];
                }
            }
        }
    }

    // Transpose.
    else if ( Teuchos::TRANS == mode )
    {
        // Make a work vector.
        Tpetra::MultiVector<double,int,SupportId> work( Y.getMap(),
                                                 Y.getNumVectors() );

        // Export X to the polynomial decomposition.
        Tpetra::Export<int,SupportId> exporter( X.getMap(), work.getMap() );
        work.doExport( X, exporter, Tpetra::INSERT );

        // Do the local mat-vec.
        Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > work_view =
            work.get2dViewNonConst();
        Teuchos::Array<double> products( poly_size * num_vec, 0.0 );
        int stride = 0;
        for ( int n = 0; n < num_vec; ++n )
        {
            stride = n*poly_size;
            for ( int p = 0; p < poly_size; ++p )
            {
                for ( int i = 0; i < local_length; ++i )
                {
                    products[stride + p] += poly_view[p][i] * work_view[n][i];
                }
            }
        }

        // Reduce the results back to the root rank.
        Teuchos::Array<double> product_sums( poly_size * num_vec, 0.0 );
#ifdef HAVE_MPI
        Teuchos::RCP< const Teuchos::MpiComm<int> > mpi_comm =
            Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( d_comm );
        Teuchos::RCP< const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm =
            mpi_comm->getRawMpiComm();
        MPI_Comm raw_comm = (*opaque_comm)();
        MPI_Reduce( products.getRawPtr(),
                    product_sums.getRawPtr(),
                    poly_size * num_vec,
                    MPI_DOUBLE,
                    MPI_SUM,
                    0,
                    raw_comm );
#else
        product_sums = products;
#endif

        // Assign the values to Y on the root rank.
        if ( 0 == d_comm->getRank() )
        {
            for ( int n = 0; n < num_vec; ++n )
            {
                for ( int p = 0; p < poly_size; ++p )
                {
                    y_view[n][p] += alpha * product_sums[n*poly_size+p];
                }
            }
        }
    }

    else
    {
        DTK_INSIST( mode == Teuchos::NO_TRANS || mode == Teuchos::TRANS );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_PolynomialMatrix.cpp
//---------------------------------------------------------------------------//

