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
 * \file   DTK_CenterDistributor.hpp
 * \author Stuart R. Slattery
 * \brief  Global center distributor.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CENTERDISTRIBUTOR_HPP
#define DTK_CENTERDISTRIBUTOR_HPP

#include "DTK_CloudDomain.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_Distributor.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class CenterDistributor
 * \brief Global source distributor.
 *
 * The CenterDistributor distributes the centers to their target
 * processes. In addition, it saves that communication plan to move source
 * field values to the same destination processes.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class CenterDistributor
{
  public:

    // Constructor.
    CenterDistributor(
	const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
	const Teuchos::ArrayView<const double>& source_centers,
	const Teuchos::ArrayView<const double>& target_centers,
	const double radius,
	Teuchos::Array<double>& target_decomp_source_centers );

    //! Destructor.
    ~CenterDistributor()
    { /* ... */ }

    // Get the number of source centers that will be distributed from this
    // process.
    int getNumExports() const
    { return d_num_exports; }

    // Get the number of source centers that will be distributed to this
    // process.
    int getNumImports() const
    { return d_num_imports; }

    // Given a set of scalar values at the given source centers in the source
    // decomposition, distribute them to the target decomposition.
    template<class T>
    void distribute( 
	const Teuchos::ArrayView<const T>& source_decomp_data,
	const Teuchos::ArrayView<T>& target_decomp_data ) const;

  private:

    // Compute the domain of the local set of centers.
    CloudDomain<DIM> localCloudDomain(
	const Teuchos::ArrayView<const double>& target_centers ) const;

  private:

    // Distributor.
    Teuchos::RCP<Tpetra::Distributor> d_distributor;

    // Number of source packets that will be imported to this process.
    int d_num_imports;

    // Number of source packets that will be exported from this process.
    int d_num_exports;

    // Expanded array of local source ids to export.
    Teuchos::Array<unsigned> d_export_ids;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_CenterDistributor_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_CENTERDISTRIBUTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_CenterDistributor.hpp
//---------------------------------------------------------------------------//

