//---------------------------------------------------------------------------//
/*!
 * \file DTK_TransferOperator_def.hpp
 * \author Stuart R. Slattery
 * \brief Transfer operator definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TRANSFEROPERATOR_DEF_HPP
#define DTK_TRANSFEROPERATOR_DEF_HPP

namespace DataTransferKit
{

namespace TransferOperator
{

/*!
 * \brief Transfer operator apply.
 * \param source_field The source field for the transfer operation.
 * \param target_field The target field for the transfer operation.
 * \param map The map for the the transfer operation.
 */
template<class SourceField, class TargetField>
void apply( const SourceField& source_field, TargetField& target_field,
	    const Teuchos::RCP<Map>& map )
{
    map->apply( source_field, target_field );
}

} // end namespace TransferOperator

} // end namespace DataTransferKit

#endif // end DTK_TRANSFEROPERATOR_DEF_HPP
