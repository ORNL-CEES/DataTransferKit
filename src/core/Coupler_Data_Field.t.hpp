//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Coupler_Data_Field.t.hpp
 * \author Stuart Slattery
 * \date   Fri Nov 18 11:57:58 2011
 * \brief  Coupler_Data_Field template member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef core_Coupler_Data_Field_t_hpp
#define core_Coupler_Data_Field_t_hpp

#include <cassert>

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR and DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * \param field_name The name of this field. Required by the
 * Coupler_Data_Source and Coupler_Data_Target interfaces to check field
 * support.
 * \param source Coupler_Data_Source implementation that will serve as the
 * data source for this field.
 * \param target Coupler_Data_Target implementation that will serve as the
 * target for this field.
 * \param scalar Set to true if this field is scalar, false if distributed.
 */
template<class DataType, class HandleType, class CoordinateType>
Data_Field<DataType,HandleType,CoordinateType>::Data_Field(
    const std::string &field_name,
    RCP_Data_Source source,
    RCP_Data_Target target,
    bool scalar)
    : d_field_name(field_name)
    , d_source(source)
    , d_target(target)
    , d_scalar(scalar)
    , d_mapped(false)
{ 
    // Require that these physics support the field being mapped.
    assert( d_source->field_supported(d_field_name) &&
	    d_target->field_supported(d_field_name) );
}

/*!
 * \brief Destructor.
 */
template<class DataType, class HandleType, class CoordinateType>
Data_Field<DataType,HandleType,CoordinateType>::~Data_Field()
{ /* ... */ }


//---------------------------------------------------------------------------//
// PUBLIC METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Assign topology maps to the field.
 * \param source_map The Tpetra map being applied to the data source for this
 * field.
 * \param target_map The Tpetra map being applied to the data target for this
 * field.
 */
template<class DataType, class HandleType, class CoordinateType>
void Data_Field<DataType,HandleType,CoordinateType>::set_mapping(
    RCP_Tpetra_Map source_map, 
    RCP_Tpetra_Map target_map)
{
    d_source_map = source_map;
    d_target_map = target_map;
    d_export = Teuchos::rcp( 
	new Tpetra::Export<HandleType>(d_source_map, d_target_map) );
    d_mapped = true;
}

} // end namespace coupler

#endif // core_Coupler_Data_Field_t_hpp

//---------------------------------------------------------------------------//
//                 end of Coupler_Data_Field.t.hpp
//---------------------------------------------------------------------------//
