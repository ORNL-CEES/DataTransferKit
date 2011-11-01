//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field/Field_DB.hh
 * \author Stuart Slattery
 * \date   Thu Oct 13 09:42:22 2011
 * \brief  Field database class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef field_Field_DB_hh
#define field_Field_DB_hh

#include "Field.hh"
#include "utils/SP.hh"
#include <string>
#include <map>

namespace field
{

//===========================================================================//
/*!
 * \class Field_DB
 * \brief The Field_DB class is a database for storing field containers in a
 * key-value pair with an associated name.
 */
/*! 
 * \example field/test/tstField_DB.cc
 *
 * Test of Field_DB.
 */
//===========================================================================//

template<class FieldType_T>
class Field_DB 
{
  public:
    
    //@{
    //! Useful typedefs.
    typedef FieldType_T                           FieldType;
    typedef Field<FieldType>                      ValueType;
    typedef denovo::SP<ValueType>                 SP_ValueType;
    typedef std::string                           KeyType;
    std::map<KeyType,SP_ValueType>                Database;
    //@}

  private:

    // Field database.
    Database d_db;

  public:
    
    //! Constructor.
    Field_DB()
    { /* ... */ }

    //! Destructor.
    ~Field_DB()
    { /* ... */ }

    //! Add a field to the database.
    inline void add_field(KeyType key);

    //! Get a field from the database.
    inline SP_ValueType get_field(KeyType key);
};

} // end namespace field

#include "Field_DB.i.hh"

#endif // field_Field_DB_hh

//---------------------------------------------------------------------------//
//              end of field/Field_DB.hh
//---------------------------------------------------------------------------//
