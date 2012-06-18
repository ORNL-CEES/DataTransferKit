//---------------------------------------------------------------------------//
/*!
 * \file DTK_FieldTools.hpp
 * \author Stuart R. Slattery
 * \brief FieldTools definition
 */
//---------------------------------------------------------------------------//

#ifndef DTK_FIELDTOOLS_HPP
#define DTK_FIELDTOOLS_HPP

#include "DTK_FieldTraits.hpp"
#include <DTK_BoundingBox.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

template<class Field>
class FieldTools
{
  public:

    //@{
    //! Typedefs. 
    typedef Field                           field_type;
    typedef FieldTraits<Field>              FT;
    typedef FT::value_type                  value_type;
    typedef Teuchos::Comm<int>              CommType;
    typedef Teuchos::RCP<const CommType>    RCP_Comm;
    //@}

    //! Constructor.
    FieldTools()
    { /* ... */ }

    //! Destructor.
    ~FieldTools()
    { /* ... */ }

    // Get the local bounding box for a coordinate field.
    BoundingBox coordLocalBoundingBox( const Field& field );

    // Get the global bounding box for a coordinate field.
    BoundingBox coordGlobalBoundingBox( const Field& field,
					const RCP_Comm& comm );

    // Get the infinity norm of a given field dimension.
    static value_type normInf( const std::size_t dim );

    // Get the L1 norm of a given field dimension.
    static value_type norm1( const std::size_t dim );

    // Get the L2 norm of a given field dimension.
    static value_type norm2( const std::size_t dim );
};

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_FieldTools_def.hpp"

#endif // end DTK_FIELDTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_FieldTools.hpp
//---------------------------------------------------------------------------//

