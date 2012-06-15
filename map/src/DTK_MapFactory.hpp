//---------------------------------------------------------------------------//
/*!
 * \file DTK_MapFactory.hpp
 * \author Stuart R. Slattery
 * \brief Map factory declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAPFACTORY_HPP
#define DTK_MAPFACTORY_HPP

#include "DTK_Map.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DataTransferKit
{

class MapFactory
{
  public:

    //@{ 
    //! Typedefs.
    typedef Teuchos::RCP<Map>    RCP_Map;
    //@}

    // Constructor.
    MapFactory();

    // Destructor.
    ~MapFactory();

    // Creation function.
    static RCP_Map create( const Teuchos::ParameterList& plist );
};

} // end namespace DataTransferKit

#endif // end DTK_MAPFACTORY_HPP

//---------------------------------------------------------------------------//
// end DTK_MapFactory.hpp
//---------------------------------------------------------------------------//


