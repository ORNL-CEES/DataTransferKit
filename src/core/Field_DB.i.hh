//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field/Field_DB.i.hh
 * \author Stuart Slattery
 * \date   Thu Oct 13 09:42:22 2011
 * \brief  Member definitions of class Field_DB.
 */
//---------------------------------------------------------------------------//
// $Id: template.i.hh,v 1.4 2008/01/04 22:50:12 9te Exp $
//---------------------------------------------------------------------------//

#ifndef field_Field_DB_i_hh
#define field_Field_DB_i_hh

namespace field
{

//---------------------------------------------------------------------------//
void Field_DB::add_field(std::string field_name)
{
    Field<ValueType> new_field;
    d_db.insert( std::make_pair(field_name,new_field) );
}

//---------------------------------------------------------------------------//
void Field_DB::get_field(std::string field_name, Field &return_field)
{
    return_field = d_db(field_name);
}

//---------------------------------------------------------------------------//

} // end namespace field

#endif // field_Field_DB_i_hh

//---------------------------------------------------------------------------//
//              end of field/Field_DB.i.hh
//---------------------------------------------------------------------------//
