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

template<class SourceMesh, class SourceDataField, class TargetDataField>
void apply( const FieldEvaluator<SourceMesh,SourceDataField>& source,
	    TargetDataField& target,
	    const Teuchos::RCP<Map>& map );

} // end namespace TransferOperator

} // end namespace DataTransferKit

#endif // DTK_TRANSFEROPERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_TransferOperator.hpp
//---------------------------------------------------------------------------//

