//---------------------------------------------------------------------------//
/*!
 * \file DTK_MapFactory.cpp
 * \author Stuart R. Slattery
 * \brief Map factory definition.
 */
//---------------------------------------------------------------------------//

#include "DTK_MapFactory.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
MapFactory::MapFactory()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
MapFactory::~MapFactory()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Creation function.
 */
MapFactory::RCP_Map MapFactory::create( const Teuchos::ParameterList& plist )
{

}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MapFactory.cpp
//---------------------------------------------------------------------------//

