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
 * \brief DTK_PointCloudOperatorFactory.hpp
 * \author Stuart R. Slattery
 * \brief Point cloud map operator factory.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POINTCLOUDOPERATORFACTORY_HPP
#define DTK_POINTCLOUDOPERATORFACTORY_HPP

#include <unordered_map>

#include "DTK_MapOperator.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class PointCloudOperatorFactory
  \brief Factory for DTK point cloud map operators.
*/
//---------------------------------------------------------------------------//
template<class Scalar>
class PointCloudOperatorFactory
{
    //! Tpetra Map typedef.
    typedef typename MapOperator<Scalar>::TpetraMap TpetraMap;

    /*!
     * \brief Constructor.
     */
    PointCloudOperatorFactory();
    
    /*!
     * \brief Creation method.
     *
     * \param domain_map Parallel map for domain vectors the created map
     * should be compatible with.
     *
     * \param range_map Parallel map for range vectors the created map should
     * be compatible with.
     *
     * \param parameters Creation parameters.
     */
    Teuchos::RCP<MapOperator<Scalar> >
    create( const Teuchos::RCP<const TpetraMap>& domain_map,
	    const Teuchos::RCP<const TpetraMap>& range_map,
	    const Teuchos::ParameterList& parameters );

  private:

    // Operator enum.
    enum PointCloudOperatorType
    {
	SPLINE_INTERPOLATION,
	MOVING_LEAST_SQUARE
    };

    // Basis enum.
    enum PointCloudBasisType
    {
	WENDLAND,
	WU,
	BUHMANN
    };

    // String name to enum map.
    std::unordered_map<std::string,int> d_name_map;

    // Basis name to enum map.
    std::unordered_map<std::string,int> d_basis_map.
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_PointCloudOperatorFactory_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_POINTCLOUDOPERATORFACTORY_HPP

//---------------------------------------------------------------------------//
// end DTK_MapOperator.hpp
//---------------------------------------------------------------------------//
