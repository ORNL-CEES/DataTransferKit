//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Data_Field.cc
 * \author Stuart Slattery
 * \date   Fri Nov 18 11:57:58 2011
 * \brief  Transfer_Data_Field member definitions.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Transfer_Data_Field.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR and DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * \param field_name The name of this field. Required by the
 * Transfer_Data_Source and Transfer_Data_Target interfaces to check field
 * support.
 * \param source Transfer_Data_Source implementation that will serve as the
 * data source for this field.
 * \param target Transfer_Data_Target implementation that will serve as the
 * target for this field.
 * \param scalar Set to true if this field is scalar, false if distributed.
 */
template<class DataType_T>
Transfer_Data_Field<DataType_T>::Transfer_Data_Field(
    const std::string &field_name,
    SP_Transfer_Data_Source source,
    SP_Transfer_Data_Target target,
    bool scalar)
    : d_field_name(field_name)
    , d_source(source)
    , d_target(target)
    , d_scalar(scalar)
    , d_mapped(false)
{ 
    // Require that these physics support the field being mapped.
    Require( d_source->field_supported(d_field_name) &&
	     d_target->field_supported(d_field_name) );
}

/*!
 * \brief Destructor.
 */
template<class DataType_T>
Transfer_Data_Field<DataType_T>::~Transfer_Data_Field()
{ /* ... */ }


//---------------------------------------------------------------------------//
// PUBLIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Assign a topology map to the field.
 * \param transfer_map The topology map being applied to this field to
 * transfer from the source to the target.
 */
template<class DataType_T>
void Transfer_Data_Field<DataType_T>::set_map(SP_Transfer_Map transfer_map)
{
    d_map = transfer_map;
    d_mapped = true;
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Transfer_Data_Field.cc
//---------------------------------------------------------------------------//
