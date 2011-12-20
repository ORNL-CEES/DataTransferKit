#ifndef wave_source_hpp
#define wave_source_hpp

#include <algorithm>

#include <Wave.hpp>

#include <Mesh_Point.hpp>
#include <Coupler_Data_Source.hpp>

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

// Data_Source interface implementation for the Wave code.
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Wave_Data_Source 
    : public Data_Source<DataType_T, HandleType_T, CoordinateType_T>
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

    // Wave object to operate on.
    RCP_Wave wave;

  public:

    Wave_Data_Source(RCP_Wave _wave)
	: wave(_wave)
    { /* ... */ }

    ~Wave_Data_Source()
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

    bool get_points(const PointType &point)
    {
	bool return_val = false;

	if ( std::find(wave->get_grid().begin(), 
		       wave->get_grid().end(),
		       point.x() )
	     != wave->get_grid().end() )
	{
	    return_val = true;
	}

	return return_val;
    }

    const Teuchos::ArrayView<DataType> send_data(const std::string &field_name)
    {
	Teuchos::ArrayView<DataType> return_view;

	if ( field_name == "WAVE_FIELD" )
	{
	    return_view = Teuchos::ArrayView<double>( wave->get_f() );
	}

	return return_view;
    }

    DataType set_global_data(const std::string &field_name)
    {
	DataType return_val = 0.0;

	return return_val;
    }
};

} // end namespace Coupler

#endif // end wave_source_hpp
