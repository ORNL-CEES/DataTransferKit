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
template<class DataType, class HandleType, class CoordinateType, int DIM>
const Teuchos::ArrayRCP<bool> 
DataSource<DataType,HandleType,CoordinateType,DIM>::are_local_points( 
    const Teuchos::ArrayView<PointType> points )
{
    Teuchos::ArrayRCP<bool> are_local( points.size(), false );
    Teuchos::ArrayRCP<bool>::iterator are_local_it;
    typename Teuchos::ArrayView<PointType>::const_iterator points_it;
    for ( points_it = points.begin(),
      are_local_it = are_local.begin(); 
	  points_it != points.end(); 
	  ++points_it, ++are_local_it )
    {
	*are_local_it = this->is_local_point( *points_it );
    }
    return are_local;
}

} // end namespace Coupler

#endif // COUPLER_DATASOURCE_DEF_HPP

//---------------------------------------------------------------------------//
//              end of Coupler_DataSource_Def.hpp
//---------------------------------------------------------------------------//
