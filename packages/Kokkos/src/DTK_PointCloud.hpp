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
 * \brief DTK_PointCloud.hpp
 * \author Stuart R. Slattery
 * \brief Point cloud.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POINTCLOUD_HPP
#define DTK_POINTCLOUD_HPP

#include <Kokkos_Core.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class PointCloud
  \brief Point cloud
 */
template<class Scalar, class GlobalOrdinal, class ExecutionSpace>
class PointCloud
{
  public:

    //@{
    //! Type aliases    
    using scalar_type = Scalar;
    using global_ordinal_type = GlobalOrdinal;
    using execution_space_type = ExecutionSpace;
    //@}

    //! Constructor.
    PointCloud( Kokkos::View<const GlobalOrdinal*,ExecutionSpace> ids,
                Kokkos::View<const Scalar**,ExecutionSpace> coordinates,
                const Teuchos::RCP<const Teuchos::Comm<int>>& comm );

    //! Get the ids.
    Kokkos::View<const GlobalOrdinal*,ExecutionSpace> ids() const;

    //! Get the coordinates.
    Kokkos::View<const Scalar**,ExecutionSpace> coordinates() const;

    //! Get the parallel communicator over which the cloud is defined.
    Teuchos::RCP<const Teuchos::Comm<int>> comm() const;

    //! Get the spatial dimension of the cloud.
    int spatialDimension() const;

    //! Get the local number of points in the cloud.
    int localSize() const;

    //! Get the global number of points in the cloud.
    GlobalOrdinal globalSize() const;
    
  private:
    
    // Global point ids. Layout = (Point)
    Kokkos::View<const GlobalOrdinal*,ExecutionSpace> _ids;

    // Point coordinates. Layout = (Point,SpaceDim)
    Kokkos::View<const Scalar**,ExecutionSpace> _coordinates;

    // Parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> _comm;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_POINTCLOUD_HPP

//---------------------------------------------------------------------------//
// end DTK_PointCloud.hpp
//---------------------------------------------------------------------------//
