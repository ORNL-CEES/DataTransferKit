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

#include "DTK_DBC.hpp"
#include "DTK_CommTools.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Ptr.hpp>

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#endif

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*
 * \brief Get MPI_COMM_WORLD in an RCP_Comm data structure if MPI is used,
 * otherwise get the default serial communicator.
 * 
 * \param comm_world World communicator.
 */
void CommTools::getCommWorld( RCP_Comm& comm_world )
{
    comm_world = Teuchos::DefaultComm<int>::getComm();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check whether two communicators own the same communication space.
 * 
 * \param comm_A Communicator A.
 *
 * \param comm_B Communicator B.
 *
 * \param comm_global An optional global communicator over which to check for
 * equality. If none is provided, MPI_COMM_WORLD will be used for an MPI
 * build. 
 *
 * \return Return true if communicator A and B operate on the same group of
 * processes. These processes do not have to have the same ranks in the
 * communicators to be considered equivalent, they are only required to be the
 * same physical process. The same result is produced independent of the
 * ordering of A and B in the input parameters.
 */
bool CommTools::equal( const RCP_Comm& comm_A, 
                       const RCP_Comm& comm_B,
                       const RCP_Comm& comm_global )
{
    RCP_Comm comm_world = comm_global;
    if ( Teuchos::is_null(comm_global) )
    {
        getCommWorld( comm_world );
    }

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
 *
 * \param comm_global An optional global communicator over which to the
 * union. If none is provided, MPI_COMM_WORLD will be used for an MPI build.
 */
void CommTools::unite( const RCP_Comm& comm_A, 
                       const RCP_Comm& comm_B,
		       RCP_Comm& comm_union,
                       const RCP_Comm& comm_global )
{
    RCP_Comm comm_world = comm_global;
    if ( Teuchos::is_null(comm_global) )
    {
        getCommWorld( comm_world );
    }

    Teuchos::Array<int> existence( comm_world->getSize(), 0 );
    Teuchos::Array<int> existence_copy( comm_world->getSize(), 0 );

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
				 &existence_copy[0] );

    int subrank;
    Teuchos::Array<int> subranks;
    Teuchos::Array<int>::const_iterator exist_begin = existence_copy.begin();
    Teuchos::Array<int>::const_iterator exist_iterator;
    for ( exist_iterator = existence_copy.begin();
	  exist_iterator != existence_copy.end();
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
 *
 * \param comm_global An optional global communicator over which to the
 * intersection. If none is provided, MPI_COMM_WORLD will be used for an MPI
 * build. 
 */
void CommTools::intersect( const RCP_Comm& comm_A, 
                           const RCP_Comm& comm_B,
			   RCP_Comm& comm_intersection,
                           const RCP_Comm& comm_global )
{
    RCP_Comm comm_world = comm_global;
    if ( Teuchos::is_null(comm_global) )
    {
        getCommWorld( comm_world );
    }

    Teuchos::Array<int> existence( comm_world->getSize(), 0 );
    Teuchos::Array<int> existence_copy( comm_world->getSize(), 0 );

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
				 &existence_copy[0] );

    int subrank;
    Teuchos::Array<int> subranks;
    Teuchos::Array<int>::const_iterator exist_begin = existence_copy.begin();
    Teuchos::Array<int>::const_iterator exist_iterator;
    for ( exist_iterator = existence_copy.begin();
	  exist_iterator != existence_copy.end();
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
/*!
 * \brief Given a comm request, check to see if it has completed.
 */
bool 
CommTools::isRequestComplete( Teuchos::RCP<Teuchos::CommRequest<int> >& handle )
{
    bool is_complete = false;

#ifdef HAVE_MPI
    DTK_REQUIRE( Teuchos::nonnull(handle) );
    Teuchos::RCP<Teuchos::MpiCommRequestBase<int> > handle_base =
	Teuchos::rcp_dynamic_cast<Teuchos::MpiCommRequestBase<int> >(handle);
    DTK_CHECK( Teuchos::nonnull(handle_base) );
    MPI_Request raw_request = handle_base->releaseRawMpiRequest();
    MPI_Status raw_status;
    int flag = 0;
    MPI_Test( &raw_request, &flag, &raw_status );
    is_complete = ( flag != 0 );
    handle = Teuchos::rcp( 
	new Teuchos::MpiCommRequestBase<int>(raw_request) );
#else
    is_complete = true;
#endif

    return is_complete;
}

//---------------------------------------------------------------------------//

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CommTools.cpp
//---------------------------------------------------------------------------//
