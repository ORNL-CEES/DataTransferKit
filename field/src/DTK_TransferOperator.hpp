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

