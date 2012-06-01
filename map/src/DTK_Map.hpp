//---------------------------------------------------------------------------//
/*!
 * \file DTK_Map.hpp
 * \author Stuart R. Slattery
 * \brief Map protocol declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAP_HPP
#define DTK_MAP_HPP

namespace DataTransferKit
{

class Map
{

  public:

    //! Constructor.
    Map()
    { /* ... */ }

    //! Destructor.
    virtual ~Map()
    { /* ... */ }

    //! Setup the map.
    virtual void setup() = 0;

    //! Apply the map.
    virtual void apply() = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_MAP_HPP

//---------------------------------------------------------------------------//
// end DTK_Map.hpp
//---------------------------------------------------------------------------//

