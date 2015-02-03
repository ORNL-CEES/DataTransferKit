//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
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
 * \brief DTK_MoabTagVector.hpp
 * \author Stuart R. Slattery
 * \brief Moab tag vector manager.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABTAGVECTOR_HPP
#define DTK_MOABTAGVECTOR_HPP

#include <vector>

#include <MBParallelComm.hpp>

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabTagVector
  \brief A class for defining Tpetra vectors over Moab tags.

  Use this class to manage Tpetra vectors and Moab tags. This class will
  create a vector over a mesh set and the tag on that mesh set. There is no
  gaurantee of consistency between the tag and the vector as the vector does
  not point to the data in the tag. Instead, the push and pull functions allow
  the user to move data between the vector and the tag as necessary.
*/
//---------------------------------------------------------------------------//
template<class Scalar>
class MoabTagVector
{
  public:

    /*!
     * \brief Constructor.
     * \param moab_mesh The mesh containing the mesh set and tag.
     * \param mesh_set The set of mesh entities over which the vector is defined.
     * \param tag The tag containing the vector data.
     */
    MoabTagVector( const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
		   const moab::EntityHandle& mesh_set,
		   const moab::Tag& tag );

    /*!
     * \brief Destructor.
     */
    ~MoabTagVector()
    { /* ... */ }

    /*!
     * \brief Get the vector over the tag.
     */
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > getVector() const
    { return d_vector; }

    /*!
     * \brief Get the vector map.
     */
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > getMap() const
    { return d_vector->getMap(); }

    /*!
     * \brief Pull data from the tag and put it into the vector.
     */
    void pullDataFromTag();

    /*!
     * \brief Push data from the vector into the tag.
     */
    void pushDataToTag();

  private:

    // The mesh over which the tag is defined.
    Teuchos::RCP<moab::ParallelComm> d_moab_mesh;

    // The mesh set over which the vector is defined.
    moab::EntityHandle d_mesh_set;

    // The tag containing the vector data.
    moab::Tag d_tag;

    // The entities over which the vector is defined.
    std::vector<moab::EntityHandle> d_entities;

    // The vector. This is a copy of the data. To put data into the vector
    // from the tag call pullDataFromTag(). To put data from the vector back
    // into the tag call pushDataToTag().
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > d_vector;  
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_MoabTagVector_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_MOABTAGVECTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabTagVector.hpp
//---------------------------------------------------------------------------//
