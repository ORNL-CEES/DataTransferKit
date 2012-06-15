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

//---------------------------------------------------------------------------//
/*!
 * \class Map
 * \brief Map protocol for data transfer.
 */
template<class SourceGeometry, class TargetGeometry,
	 class SourceField, class TargetField>
class Map
{
  public:

    //! Constructor.
    Map()
    { /* ... */ }

    //! Destructor.
    virtual ~Map()
    { /* ... */ }

    //! Setup the map for a given geometry.
    virtual void setup( const SourceGeometry& source_geometry,
			const TargetGeometry& target_geometry ) = 0;

    //! Apply the map for a given field.
    virtual void apply( const SourceField& source_field,
			TargetField& target_field ) = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_MAP_HPP

//---------------------------------------------------------------------------//
// end DTK_Map.hpp
//---------------------------------------------------------------------------//

