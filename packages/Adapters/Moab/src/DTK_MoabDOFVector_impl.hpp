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
 * \brief DTK_MoabDOFVector_impl.hpp
 * \author Stuart R. Slattery
 * \brief Moab mesh tag DOF vector.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABDOFVECTOR_IMPL_HPP
#define DTK_MOABDOFVECTOR_IMPL_HPP

#include <vector>

#include "DTK_DBC.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Map.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Given a Moab tag, create a Tpetra vector that maps to the tag DOFs on the
// given mesh set and pull the data from the tag.
template<class Scalar>
Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > 
MoabDOFVector::pullTpetraMultiVectorFromMoabTag( 
    const moab::ParallelComm& moab_mesh,
    const moab::EntityHandle& mesh_set,
    const moab::Tag& tag )
{
    // Get the dimension of the tag.
    int tag_dim = 0;
    DTK_CHECK_ERROR_CODE(
	moab_mesh.get_moab()->tag_get_length( tag, tag_dim )
	);

    // Get the entities in the set.
    std::vector<moab::EntityHandle> entities;
    DTK_CHECK_ERROR_CODE(
	moab_mesh.get_moab()->get_entities_by_handle( mesh_set, entities )
	);
    int num_entities = entities.size();

    // Extract the MPI communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::rcp( new Teuchos::MpiComm<int>(moab_mesh.comm()) );

    // Build a map. The DOF ids are the entity handles.
    Teuchos::Array<std::size_t> dof_ids( entities.begin(), entities.end() );
    Teuchos::RCP<const Tpetra::Map<int,std::size_t> > map =
	Tpetra::createNonContigMap<int,std::size_t>( dof_ids(), comm );

    // Create a Tpetra vector.
    Teuchos::RCP<Tpetra::MultiVector<Scalar,int,std::size_t> > vector =
	Tpetra::createMultiVector<Scalar,int,std::size_t>( map, tag_dim );

    // Only populate the vector if there are entities locally in the set.
    if ( 0 < num_entities )
    {
	// Extract pointers to the data for each entity handle.
	Teuchos::Array<const void*> data_ptrs( num_entities );
	DTK_CHECK_ERROR_CODE(
	    moab_mesh.get_moab()->tag_get_by_ptr( 
		tag,
		entities.data(),
		num_entities,
		data_ptrs.getRawPtr() )
	    );

	// Extract the data in blocked form.
	Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > vector_view =
	    vector->get2dViewNonConst();
	const Scalar* entity_data;
	for ( int e = 0; e < num_entities; ++e )
	{
	    entity_data = static_cast<const Scalar*>( data_ptrs[e] );
	    for ( int d = 0; d < tag_dim; ++d )
	    {
		vector_view[d][e] = entity_data[d];
	    }
	}
    }

    // Return the vector.
    return vector;    
}

//---------------------------------------------------------------------------//
// Given a Tpetra vector of DOF data, push the data into a given Moab tag.
template<class Scalar>
void MoabDOFVector::pushTpetraMultiVectorToMoabTag(
    const Tpetra::MultiVector<Scalar,int,std::size_t>& tag_dofs,
    const moab::ParallelComm& moab_mesh,
    const moab::EntityHandle& mesh_set,
    const moab::Tag& tag )
{
    // Get the entities in the set.
    std::vector<moab::EntityHandle> entities;
    DTK_CHECK_ERROR_CODE(
    	moab_mesh.get_moab()->get_entities_by_handle( mesh_set, entities )
    	);
    int num_entities = entities.size();

    if ( 0 < num_entities )
    {
	// Get the dimension of the tag.
	int tag_dim = 0;
	DTK_CHECK_ERROR_CODE(
	    moab_mesh.get_moab()->tag_get_length( tag, tag_dim )
	    );
	DTK_CHECK( tag_dim == Teuchos::as<int>(tag_dofs.getNumVectors()) );

	// Extract data from the vector in interleaved format.
	Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > vector_view =
	    tag_dofs.get2dView();
	Teuchos::Array<Scalar> interleaved_data( num_entities * tag_dim );
	for ( int d = 0; d < tag_dim; ++d )
	{
	    for ( int e = 0; e < num_entities; ++e )
	    {
		interleaved_data[ e*tag_dim + d ] = vector_view[d][e];
	    }
	}

	// Set the data with the tag.
	DTK_CHECK_ERROR_CODE(
	    moab_mesh.get_moab()->tag_set_data( 
		tag,
		entities.data(),
		num_entities,
		static_cast<void*>(interleaved_data.getRawPtr()) )
	    );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABDOFVECTOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabDOFVector_impl.hpp
//---------------------------------------------------------------------------//
