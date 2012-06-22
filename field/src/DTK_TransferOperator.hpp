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
 * \file DTK_TransferOperator.hpp
 * \author Stuart R. Slattery
 * \brief Transfer operator declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TRANSFEROPERATOR_HPP
#define DTK_TRANSFEROPERATOR_HPP

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

template<class Map>
class TransferOperator
{
  public:

    //@{
    //! Typedefs.
    typedef Map                  map_type;
    typedef Teuchos::RCP<Map>    RCP_Map;

    //! Constructor.
    TransferOperator( const RCP_Map& map )
	: d_map( map )
    { /* ... */ }

    //! Destructor.
    ~TransferOperator()
    { /* ... */ }

    // Transfer operator setup.
    template<class SourceGeometry, class TargetGeometry>
    inline void setup( const SourceGeometry& source_geometry,
		       const TargetGeometry& target_geometry );

    // Transfer operator apply.
    template<class SourceEvaluator, class TargetSpace>
    inline void apply( const SourceEvaluator& source_evaluator, 
		       TargetSpace& target_space );

  private:

    // Map.
    RCP_Map d_map;
};

//---------------------------------------------------------------------------//
// Inline functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Transfer operator setup.
 * \param source_geometry The source geometry for the transfer operation.
 * \param target_geometry The target geometry for the transfer operation.
 */
template<class Map>
template<class SourceGeometry, class TargetGeometry>
void TransferOperator<Map>::setup( const SourceGeometry& source_geometry, 
				   const TargetGeometry& target_geometry )
{
    d_map->setup( source_geometry, target_geometry );
}
 
//---------------------------------------------------------------------------//
/*!
 * \brief Transfer operator apply.
 * \param source_evaluator The source field evaluator for the transfer
 * operation. 
 * \param target_space The target field space for the transfer operation.
 */
template<class Map>
template<class SourceEvaluator, class TargetSpace>
void TransferOperator<Map>::apply( const SourceEvaluator& source_evaluator, 
				   TargetSpace& target_space )
{
    d_map->apply( source_evaluator, target_space );
}
 
//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // DTK_TRANSFEROPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_TransferOperator.hpp
//---------------------------------------------------------------------------//

