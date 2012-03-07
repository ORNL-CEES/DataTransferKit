//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Coupler_DataSource_Def.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:43 2011
 * \brief  Default implementation for data source applications.
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_DATASOURCE_DEF_HPP
#define COUPLER_DATASOURCE_DEF_HPP

namespace Coupler
{

// Default implementation.
template<class DataType, class HandleType, class CoordinateType>
const Teuchos::ArrayRCP<bool> 
DataSource<DataType,HandleType,CoordinateType>::are_local_points( 
    const Teuchos::ArrayView<PointType> points )
{
    Teuchos::ArrayRCP<bool> are_local( points.size(), false );
    Teuchos::ArrayRCP<bool>::iterator are_local_it = are_local.begin();
    typename Teuchos::ArrayView<PointType>::const_iterator point_it;
    for ( point_it = points.begin(); 
	  point_it != points.end(); 
	  ++point_it, ++are_local_it )
    {
	*are_local_it = this->is_local_point( *point_it );
    }
    return are_local;
}

} // end namespace Coupler

#endif // COUPLER_DATASOURCE_DEF_HPP

//---------------------------------------------------------------------------//
//              end of Coupler_DataSource_Def.hpp
//---------------------------------------------------------------------------//
