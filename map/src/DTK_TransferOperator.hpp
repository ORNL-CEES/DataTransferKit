//---------------------------------------------------------------------------//
/*!
 * \file DTK_TransferOperator.hpp
 * \author Stuart R. Slattery
 * \brief Transfer operator declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TRANSFEROPERATOR_HPP
#define DTK_TRANSFEROPERATOR_HPP

#include <DTK_Map.hpp>

namespace DataTransferKit
{

namespace TransferOperator
{

// Transfer operator apply.
template<class SourceField, class TargetField>
void apply( const SourceField& source_field, TargetField& target_field,
	    const Teuchos::RCP<Map>& map );

} // end namespace TransferOperator

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_TransferOperator_def.hpp"

#endif // DTK_TRANSFEROPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_TransferOperator.hpp
//---------------------------------------------------------------------------//

