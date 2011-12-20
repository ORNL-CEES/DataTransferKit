#ifndef wave_target_hpp
#define wave_target_hpp

#include <algorithm>

#include "Wave.hpp"

#include <Mesh_Point.hpp>
#include <Coupler_Data_Target.hpp>

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"

//---------------------------------------------------------------------------//
// Get the current default communicator.
template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
namespace Coupler {

// Data_Target interface implementation for the Wave code.
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Wave_Data_Target 
    : public Data_Target<DataType_T, HandleType_T, CoordinateType_T>
{
  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef Point<HandleType,CoordinateType>         PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    typedef Teuchos::RCP<Wave>                       RCP_Wave;

  private:

    RCP_Wave wave;

  public:

    Wave_Data_Target(RCP_Wave _wave)
	: wave(_wave)
    { /* ... */ }

    ~Wave_Data_Target()
    { /* ... */ }

    RCP_Communicator comm()
    {
	return getDefaultComm<OrdinalType>();
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
	    std::vector<PointType> local_points;
	    Teuchos::ArrayView<const double> local_grid( wave->get_grid() );
	    Teuchos::ArrayView<double>::const_iterator grid_it;
	    int n = 0;
	    int global_handle;
	    for (grid_it = local_grid.begin(); 
		 grid_it != local_grid.end();
		 ++grid_it, ++n)
	    {
		global_handle = getDefaultComm<int>()->getSize()*
				getDefaultComm<int>()->getRank() + n;
		local_points.push_back( 
		    PointType( handle, *grid_it, 0.0, 0.0) );
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
	    return_view = Teuchos::ArrayView<DataType>( wave->set_damping() );
	}

	return return_view;
    }

    void get_global_data(const std::string &field_name,
			 const DataType &data)
    { /* ... */ }
};

//---------------------------------------------------------------------------//

} // end namespace Coupler

#endif // end wave_target_hpp
