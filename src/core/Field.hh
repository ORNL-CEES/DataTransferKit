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
    //@}

  private:

    

  public:

    //! Constructor.
    Field()
    { /* ... */ }

    //! Destructor.
    ~Field()
    { /* ... */ }

};

#include "Field.i.hh"

} // end namespace field

#endif // field_Field_hh

//---------------------------------------------------------------------------//
//              end of field/Field.hh
//---------------------------------------------------------------------------//
