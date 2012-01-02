#ifndef damper_target_hpp
#define damper_target_hpp

#include <algorithm>
#include <vector>

#include "Damper.hpp"

#include <Mesh_Point.hpp>
#include <Coupler_DataTarget.hpp>

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"

//---------------------------------------------------------------------------//
namespace Coupler {

// DataTarget interface implementation for the Damper code.
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Damper_DataTarget 
    : public DataTarget<DataType_T, HandleType_T, CoordinateType_T>
{
  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef Point<HandleType,CoordinateType>         PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    typedef Teuchos::RCP<Damper>                       RCP_Damper;

  private:

    RCP_Damper damper;
    std::vector<PointType> local_points;

  public:

    Damper_DataTarget(RCP_Damper _damper)
	: damper(_damper)
    { /* ... */ }

    ~Damper_DataTarget()
    { /* ... */ }

    RCP_Communicator comm()
    {
	return damper->get_comm();
    }

    bool field_supported(const std::string &field_name)
    {
	bool return_val = false;

	if (field_name == "WAVE_FIELD")
	{
	    return_val = true;
	}

	return return_val;
    }

    const Teuchos::ArrayView<PointType> 
    set_points(const std::string &field_name)
    {
	Teuchos::ArrayView<PointType> return_view;

	if ( field_name == "WAVE_FIELD" )
	{
	    local_points.clear();
	    Teuchos::ArrayView<const DataType> local_grid( damper->get_grid() );
	    Teuchos::ArrayView<DataType>::const_iterator grid_it;
	    int n = 0;
	    int global_handle;
	    for (grid_it = local_grid.begin(); 
		 grid_it != local_grid.end();
		 ++grid_it, ++n)
	    {
		global_handle = damper->get_comm()->getRank() *
				local_grid.size() + n;
		local_points.push_back( 
		    PointType( global_handle, *grid_it, 0.0, 0.0) );
	    }
	    return_view = Teuchos::ArrayView<PointType>(local_points);
	}

	return return_view;
    }

    Teuchos::ArrayView<DataType> receive_data(const std::string &field_name)
    {
	Teuchos::ArrayView<DataType> return_view;

	if ( field_name == "WAVE_FIELD" )
	{
	    return_view = Teuchos::ArrayView<DataType>( damper->set_wave_data() );
	}

	return return_view;
    }

    void get_global_data(const std::string &field_name,
			 const DataType &data)
    { /* ... */ }
};

//---------------------------------------------------------------------------//

} // end namespace Coupler

#endif // end damper_target_hpp
