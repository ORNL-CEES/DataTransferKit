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
 * \file DTK_CommTools.cpp
 * \author Stuart R. Slattery
 * \brief CommTools definition.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "DTK_CommTools.hpp"

#include <mpi.h>

#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*
 * \brief Get MPI_COMM_WORLD in an RCP_Comm data structure.
 * 
 * \param mpi_comm_world MPI_COMM_WORLD wrapped in Teuchos::Comm object.
 */
void CommTools::getMpiCommWorld( RCP_Comm& mpi_comm_world )
{
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm =
	Teuchos::rcp( new Teuchos::OpaqueWrapper<MPI_Comm>( MPI_COMM_WORLD ) );
    mpi_comm_world = Teuchos::rcp( new Teuchos::MpiComm<int>( opaque_comm ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check whether two communicators own the same communication space.
 * 
 * \param comm_A Communicator A.
 *
 * \param comm_B Communicator B.
 *
 * \return Return true if communicator A and B operate on the same group of
 * processes. These processes do not have to have the same ranks in the
 * communicators to be considered equivalent, they are only required to be the
 * same physical process. The same result is produced independent of the
 * ordering of A and B in the input parameters.
 */
bool CommTools::equal( const RCP_Comm& comm_A, const RCP_Comm& comm_B )
{
    RCP_Comm comm_world;
    getMpiCommWorld( comm_world );

    int existence = 0;

    if ( !comm_A.is_null() )
    {
	++existence;
    }

    if ( !comm_B.is_null() )
    {
	++existence;
    }

    int local_not_equal = 0;
    if ( existence == 1 )
    {
	local_not_equal = 1;
    }

    int global_not_equal = 0;
    Teuchos::reduceAll<int,int>( *comm_world, 
				 Teuchos::REDUCE_SUM,
				 local_not_equal, 
				 Teuchos::Ptr<int>(&global_not_equal) );
    
    if ( global_not_equal > 0 )
    {
	return false;
    }
    
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the union of two communicators.
 * 
 * \param comm_A Communicator A.
 *
 * \param comm_B Communicator B.
 *
 * \param comm_union Return the union of communicator A and B. The union of A
 * and B is defined as the physical processes that exist in either
 * communicator. The same result is produced independent of the ordering of A
 * and B in the input parameters. 
 */
void CommTools::unite( const RCP_Comm& comm_A, const RCP_Comm& comm_B,
		       RCP_Comm& comm_union )
{
    RCP_Comm comm_world;
    getMpiCommWorld( comm_world );

    Teuchos::Array<int> existence( comm_world->getSize(), 0 );

    if ( !comm_A.is_null() )
    {
	existence[ comm_world->getRank() ] += 1;
    }

    if ( !comm_B.is_null() )
    {
	existence[ comm_world->getRank() ] += 1;
    }
    comm_world->barrier();

    Teuchos::reduceAll<int,int>( *comm_world,
				 Teuchos::REDUCE_SUM,
				 (int) existence.size(),
				 &existence[0],
				 &existence[0] );

    int subrank;
    Teuchos::Array<int> subranks;
    Teuchos::Array<int>::const_iterator exist_begin = existence.begin();
    Teuchos::Array<int>::const_iterator exist_iterator;
    for ( exist_iterator = existence.begin();
	  exist_iterator != existence.end();
	  ++exist_iterator )
    {
	if ( *exist_iterator > 0 )
	{
	    subrank = std::distance( exist_begin, exist_iterator );
	    subranks.push_back( subrank );
	}
    }
   
    comm_union = comm_world->createSubcommunicator( subranks() );    
}

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the intersection of two communicators.
 * 
 * \param comm_A Communicator A.
 *
 * \param comm_B Communicator B.
 *
 * \param comm_intersection Return the intersection of communicator A and
 * B. The intersection of A and B is defined to be the group of physical
 * processes that exist in both A and B. The same result is produced
 * independent of the ordering of A and B in the input parameters.
 */
void CommTools::intersect( const RCP_Comm& comm_A, const RCP_Comm& comm_B,
			   RCP_Comm& comm_intersection )
{
    RCP_Comm comm_world;
    getMpiCommWorld( comm_world );

    Teuchos::Array<int> existence( comm_world->getSize(), 0 );

    if ( !comm_A.is_null() )
    {
	existence[ comm_world->getRank() ] += 1;
    }

    if ( !comm_B.is_null() )
    {
	existence[ comm_world->getRank() ] += 1;
    }

    Teuchos::reduceAll<int,int>( *comm_world,
				 Teuchos::REDUCE_SUM,
				 (int) existence.size(),
				 &existence[0],
				 &existence[0] );

    int subrank;
    Teuchos::Array<int> subranks;
    Teuchos::Array<int>::const_iterator exist_begin = existence.begin();
    Teuchos::Array<int>::const_iterator exist_iterator;
    for ( exist_iterator = existence.begin();
	  exist_iterator != existence.end();
	  ++exist_iterator )
    {
	if ( *exist_iterator == 2)
	{
	    subrank = std::distance( exist_begin, exist_iterator );
	    subranks.push_back( subrank );
	}
    }
   
    comm_intersection = comm_world->createSubcommunicator( subranks() );    
}

//---------------------------------------------------------------------------//

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CommTools.cpp
//---------------------------------------------------------------------------//
