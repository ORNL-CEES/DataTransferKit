//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field/Field.i.hh
 * \author Stuart Slattery
 * \date   Thu Oct 13 09:42:01 2011
 * \brief  Member definitions of class Field.
 */
//---------------------------------------------------------------------------//
// $Id: template.i.hh,v 1.4 2008/01/04 22:50:12 9te Exp $
//---------------------------------------------------------------------------//

#ifndef field_Field_i_hh
#define field_Field_i_hh

namespace field
{

//---------------------------------------------------------------------------//
// Set the domain of the field.
void Field::set_domain(Const_Iterator begin, Const_Iterator end)
{
    d_domain_begin = begin;
    d_domain_end = end;
}

//---------------------------------------------------------------------------//
// Set the range of the field.
void Field::set_range(Iterator begin, Iterator end)
{
    d_range_begin = begin;
    d_range_end = end;
}

//---------------------------------------------------------------------------//
// Get the domain of the field.
void Field::get_domain(Const_Iterator &begin, Const_Iterator &end)
{
    begin = d_domain_begin;
    end = d_domain_end;
}

//---------------------------------------------------------------------------//
// Set the range of the field.
void Field::set_range(Iterator &begin, Iterator &end)
{
    begin = d_range_begin;
    end = d_range_end;
}

//---------------------------------------------------------------------------//

} // end namespace field

#endif // field_Field_i_hh

//---------------------------------------------------------------------------//
//              end of field/Field.i.hh
//---------------------------------------------------------------------------//
