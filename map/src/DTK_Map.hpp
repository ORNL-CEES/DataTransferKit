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
class Map
{

  public:

    //@{
    //! Typedefs.
    //@}

    //! Constructor.
    Map()
    { /* ... */ }

    //! Destructor.
    virtual ~Map()
    { /* ... */ }

    //! Setup the map for a given geometry.
    template<class SourceGeometry, class TargetGeometry>
    virtual void setup( const SourceGeometry& source_geom,
			const TargetGeometry& target_geom ) = 0;

    //! Apply the map for a given field.
    template<class SourceField, class TargetField>
    virtual void apply( const SourceField& source_field,
			TargetField& target_field ) = 0;
};

} // end namespace DataTransferKit

#endif // end DTK_MAP_HPP

//---------------------------------------------------------------------------//
// end DTK_Map.hpp
//---------------------------------------------------------------------------//

