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

#include <vector>

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
    typedef typename FieldType::value_type         ValueType;
    typedef ValueType::iterator                    Iterator;
    typedef ValueType::const_iterator              Const_Iterator;
    //@}

  private:
    
    //@{
    //! Domain iterators.
    Const_Iterator d_domain_begin;
    Const_Iterator d_domain_end;
    //@}
    
    //@{
    //! Range iterators.
    Iterator d_range_begin;
    Iterator d_range_end;
    //@}

  public:

    //! Constructor.
    Field()
    { /* ... */ }

    //! Destructor.
    ~Field()
    { /* ... */ }

    //! Set the domain of the field.
    inline void set_domain(Const_Iterator begin, Const_Iterator end);

    //! Set the range of the field.
    inline void set_range(Iterator begin, Iterator end);

    //! Get the domain of the field.
    inline void get_domain(Const_Iterator &begin, Const_Iterator &end);

    //! Set the range of the field.
    inline void set_range(Iterator &begin, Iterator &end);    
};

#include "Field.i.hh"

} // end namespace field

#endif // field_Field_hh

//---------------------------------------------------------------------------//
//              end of field/Field.hh
//---------------------------------------------------------------------------//
