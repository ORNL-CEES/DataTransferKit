//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field/Field.hh
 * \author Stuart Slattery
 * \date   Thu Oct 13 09:42:01 2011
 * \brief  Field class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef field_Field_hh
#define field_Field_hh

#include "Transfer_Map.hh"

namespace field
{

//===========================================================================//
/*!
 * \class Field
 * \brief The field class is a container for data fields.
 */
/*! 
 * \example field/test/tstField.cc
 *
 * Test of Field.
 */
//===========================================================================//

template<class FieldType_T>
class Field 
{
  public:
    
    //@{
    //! Useful typedefs.
    typedef FieldType_T                            FieldType;
    typedef int                                    Handle;
    typedef int                                    Ordinate;
    typedef Transfer_Map<Handle,Ordinate>          Map;
    typedef nemesis::SP<Map>                       SP_Transfer_Map;
    //@}

  private:

    // Topology map for transfer from A to B.
    SP_Transfer_Map d_map_a2b;

    // Topology map for transfer from B to A.
    SP_Transfer_Map d_map_b2a;

  public:

    //! Constructor.
    Field();

    //! Destructor.
    ~Field();

    //! Get the map for A to B.
    SP_Transfer_Map map_a2b();

    //! Get the map for B to A.
    SP_Transfer_Map map_b2a();
};

#include "Field.i.hh"

} // end namespace field

#endif // field_Field_hh

//---------------------------------------------------------------------------//
//              end of field/Field.hh
//---------------------------------------------------------------------------//
