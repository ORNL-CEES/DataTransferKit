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
 * \file tstPOD_PointCloudEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief POD_PointCloudEntity unit tests.
 */
//---------------------------------------------------------------------------//
#include "DTK_POD_PointCloudEntityIterator.hpp"
#include "DTK_POD_Types.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( POD_PointCloudEntityIterator, blocked_test )
{
  // get the raw mpi communicator
  Teuchos::RCP<const Teuchos::Comm<int>> teuchos_comm =
    Teuchos::DefaultComm<int>::getComm();
  MPI_Comm mpi_comm =
    *Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(
      teuchos_comm )->getRawMpiComm();
  int const comm_rank = teuchos_comm->getRank();

  // Build a blocked point cloud.
  std::srand(123*comm_rank);
  int const space_dim = 3;
  unsigned const num = 3000;
  Teuchos::Array<double> coord(space_dim*num);
  Teuchos::Array<DataTransferKit::EntityId> global_ids(num);
  for ( unsigned i = 0; i < num; ++i )
  {
      coord[i+0*num] = (double) std::rand() / (double) RAND_MAX + comm_rank;
      coord[i+1*num] = (double) std::rand() / (double) RAND_MAX;
      coord[i+2*num] = (double) std::rand() / (double) RAND_MAX;

      global_ids[i] = num*comm_rank + i;
  }

  // Build an iterator.
  auto select_all = [](DataTransferKit::Entity e){ return true; };
  DataTransferKit::EntityIterator iterator =
      DataTransferKit::POD_PointCloudEntityIterator( coord.getRawPtr(),
                                                     global_ids.getRawPtr(),
                                                     num,
                                                     space_dim,
                                                     DataTransferKit::BLOCKED,
                                                     comm_rank,
                                                     select_all );

  // Check the iterator.
  TEST_EQUALITY( num, iterator.size() );

  // Loop through the point cloud and check entities.
  auto begin_it = iterator.begin();
  auto end_it = iterator.end();
  unsigned i = 0;
  for ( iterator = begin_it; iterator != end_it; ++iterator, ++i )
  {
      TEST_EQUALITY( iterator->id(), global_ids[i] );
      TEST_EQUALITY( iterator->ownerRank(), comm_rank );
      TEST_EQUALITY( iterator->topologicalDimension(), 0 );
      TEST_EQUALITY( iterator->physicalDimension(), space_dim );

      Teuchos::Tuple<double,6> box;
      iterator->boundingBox( box );
      TEST_EQUALITY( box[0], coord[i] );
      TEST_EQUALITY( box[1], coord[i + 1*num] );
      TEST_EQUALITY( box[2], coord[i + 2*num] );
      TEST_EQUALITY( box[3], coord[i] );
      TEST_EQUALITY( box[4], coord[i + 1*num] );
      TEST_EQUALITY( box[5], coord[i + 2*num] );

      TEST_ASSERT( !iterator->inBlock(1) );
      TEST_ASSERT( !iterator->onBoundary(1) );
  }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( POD_PointCloudEntityIterator, interleaved_test )
{
  // get the raw mpi communicator
  Teuchos::RCP<const Teuchos::Comm<int>> teuchos_comm =
    Teuchos::DefaultComm<int>::getComm();
  MPI_Comm mpi_comm =
    *Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(
      teuchos_comm )->getRawMpiComm();
  int const comm_rank = teuchos_comm->getRank();

  // Build a interleaved point cloud.
  std::srand(123*comm_rank);
  int const space_dim = 3;
  unsigned const num = 3000;
  Teuchos::Array<double> coord(space_dim*num);
  Teuchos::Array<DataTransferKit::EntityId> global_ids(num);
  for ( unsigned i = 0; i < num; ++i )
  {
      coord[space_dim*i + 0] = (double) std::rand() / (double) RAND_MAX + comm_rank;
      coord[space_dim*i + 1] = (double) std::rand() / (double) RAND_MAX;
      coord[space_dim*i + 2] = (double) std::rand() / (double) RAND_MAX;

      global_ids[i] = num*comm_rank + i;
  }

  // Build an iterator.
  auto select_all = [](DataTransferKit::Entity e){ return true; };
  DataTransferKit::EntityIterator iterator =
      DataTransferKit::POD_PointCloudEntityIterator( coord.getRawPtr(),
                                                     global_ids.getRawPtr(),
                                                     num,
                                                     space_dim,
                                                     DataTransferKit::INTERLEAVED,
                                                     comm_rank,
                                                     select_all );

  // Check the iterator.
  TEST_EQUALITY( num, iterator.size() );

  // Loop through the point cloud and check entities.
  auto begin_it = iterator.begin();
  auto end_it = iterator.end();
  unsigned i = 0;
  for ( iterator = begin_it; iterator != end_it; ++iterator, ++i )
  {
      TEST_EQUALITY( iterator->id(), global_ids[i] );
      TEST_EQUALITY( iterator->ownerRank(), comm_rank );
      TEST_EQUALITY( iterator->topologicalDimension(), 0 );
      TEST_EQUALITY( iterator->physicalDimension(), space_dim );

      Teuchos::Tuple<double,6> box;
      iterator->boundingBox( box );
      TEST_EQUALITY( box[0], coord[space_dim*i + 0] );
      TEST_EQUALITY( box[1], coord[space_dim*i + 1] );
      TEST_EQUALITY( box[2], coord[space_dim*i + 2] );
      TEST_EQUALITY( box[3], coord[space_dim*i + 0] );
      TEST_EQUALITY( box[4], coord[space_dim*i + 1] );
      TEST_EQUALITY( box[5], coord[space_dim*i + 2] );

      TEST_ASSERT( !iterator->inBlock(1) );
      TEST_ASSERT( !iterator->onBoundary(1) );
  }
}

//---------------------------------------------------------------------------//
// end tstPOD_PointCloudEntityIterator.cpp
//---------------------------------------------------------------------------//

