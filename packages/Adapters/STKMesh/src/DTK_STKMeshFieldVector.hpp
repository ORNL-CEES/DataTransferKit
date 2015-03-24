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
 * \brief DTK_STKMeshFieldVector.hpp
 * \author Stuart R. Slattery
 * \brief STK field vector manager.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHFIELDVECTOR_HPP
#define DTK_STKMESHFIELDVECTOR_HPP

#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshFieldVector
  \brief A class for defining Tpetra vectors over STK fields.

  Use this class to manage Tpetra vectors and STK fields. This class will
  create a vector over a mesh set and the field on that mesh set. There is no
  gaurantee of consistency between the field and the vector as the vector does
  not point to the data in the field. Instead, the push and pull functions allow
  the user to move data between the vector and the field as necessary.
*/
//---------------------------------------------------------------------------//
template<class Scalar,class FieldType>
class STKMeshFieldVector
{
  public:

    /*!
     * \brief Constructor.
     */
    STKMeshFieldVector( const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
			const Teuchos::Ptr<FieldType>& field,
			const int field_dim );

    /*!
     * \brief Get the vector over the field.
     */
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > getVector() const
    { return d_vector; }

    /*!
     * \brief Get the vector map.
     */
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > getMap() const
    { return d_vector->getMap(); }

    /*!
     * \brief Pull data from the field and put it into the vector.
     */
    void pullDataFromField();

    /*!
     * \brief Push data from the vector into the field.
     */
    void pushDataToField();

  private:

    // The mesh over which the field is defined.
    Teuchos::RCP<stk::mesh::BulkData> d_bulk_data;

    // The field containing the vector data.
    Teuchos::Ptr<FieldType> d_field;

    // The dimension of the field.
    int d_field_dim;

    // The enitities over which the field is defined.
    std::vector<stk::mesh::Entity> d_field_entities;

    // The vector. This is a copy of the data. To put data into the vector
    // from the field call pullDataFromField(). To put data from the vector back
    // into the field call pushDataToField().
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > d_vector;  
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_STKMeshFieldVector_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_STKMESHFIELDVECTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshFieldVector.hpp
//---------------------------------------------------------------------------//
