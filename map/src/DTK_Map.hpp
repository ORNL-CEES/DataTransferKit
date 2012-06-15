//---------------------------------------------------------------------------//
/*!
 * \file DTK_Map.hpp
 * \author Stuart R. Slattery
 * \brief Map protocol declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAP_HPP
#define DTK_MAP_HPP

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class Map
 * \brief Map protocol for data transfer.
 */
template<typename GlobalOrdinal>
class Map
{

  public:

    //@{
    //! Typedefs.
    typedef GlobalOrdinal                   global_ordinal_type;
    typedef Tpetra::Map<GlobalOrdinal>      TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>   RCP_TpetraMap;
    //@}

    //! Constructor.
    Map()
    { /* ... */ }

    //! Destructor.
    virtual ~Map()
    { /* ... */ }

    //! Build the map.
    virtual void build( RCP_TpetraMap& export_map,
			RCP_TpetraMap& import_map ) = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_MAP_HPP

//---------------------------------------------------------------------------//
// end DTK_Map.hpp
//---------------------------------------------------------------------------//

