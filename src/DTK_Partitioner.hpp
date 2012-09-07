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
 * \file DTK_Partitioner.hpp
 * \author Stuart R. Slattery
 * \brief Partitioner interface definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PARTITIONER_HPP
#define DTK_PARTITIONER_HPP

#include "DTK_BoundingBox.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class Partitioner
 * \brief Partitioner base class.
 */
//---------------------------------------------------------------------------//
class Partitioner
{
  public:
    
    // Constructor.
    Partitioner()
    { /* ... */ }

    // Destructor.
    virtual ~Partitioner()
    { /* ... */ }

    // Compute partitioning of the mesh.
    virtual void partition() = 0;

    // Get the destination process for a point given its coordinates.
    virtual int 
    getPointDestinationProc( Teuchos::Array<double> coords ) const = 0;

    // Get the destination processes for a bounding box.
    virtual Teuchos::Array<int> 
    getBoxDestinationProcs( const BoundingBox& box ) const = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_PARTITIONER_HPP

//---------------------------------------------------------------------------//
// end DTK_Partitioner.hpp
//---------------------------------------------------------------------------//
