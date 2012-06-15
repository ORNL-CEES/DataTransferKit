//---------------------------------------------------------------------------//
/*!
 * \file DTK_TransferOperator.hpp
 * \author Stuart R. Slattery
 * \brief Transfer operator declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TRANSFEROPERATOR_HPP
#define DTK_TRANSFEROPERATOR_HPP

namespace DataTransferKit
{

namespace TransferOperator
{

//! Transfer operator apply.
template<class SourceMesh, class SourceDataField, 
	 class TargetDataField, class GlobalOrdinal>
void apply( const FieldEvaluator<SourceMesh,SourceDataField>& source,
	    TargetDataField& target,
	    const Teuchos::RCP< Map<GlobalOrdinal> >& map );

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

