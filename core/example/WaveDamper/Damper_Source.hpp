#ifndef damper_source_hpp
#define damper_source_hpp

#include <algorithm>

#include "Damper.hpp"

#include <Coupler_Point.hpp>
#include <Coupler_DataSource.hpp>

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"

//---------------------------------------------------------------------------//
namespace Coupler {

// DataSource interface implementation for the Damper code.
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Damper_DataSource 
    : public DataSource<DataType_T, HandleType_T, CoordinateType_T>
{
  public:

    typedef double                                   DataType;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    typedef int                                      OrdinalType;
    typedef Point<HandleType,CoordinateType>         PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    typedef Teuchos::RCP<Damper>                     RCP_Damper;

  private:

    // Damper object to operate on.
    RCP_Damper damper;

  public:

    Damper_DataSource(RCP_Damper _damper)
	: damper(_damper)
    { /* ... */ }

    ~Damper_DataSource()
    { /* ... */ }

    RCP_Communicator get_source_comm()
    {
	return damper->get_comm();
    }

    bool is_field_supported(const std::string &field_name)
    {
	bool return_val = false;

	if (field_name == "DAMPER_SOURCE_FIELD")
	{
	    return_val = true;
	}

	return return_val;
    }

    bool is_local_point(const PointType &point)
    {
	bool return_val = false;

	if ( std::find(damper->get_grid().begin(), 
		       damper->get_grid().end(),
		       point.x() )
	     != damper->get_grid().end() )
	{
	    return_val = true;
	}

	return return_val;
    }

    const Teuchos::ArrayView<DataType> 
    get_source_data(const std::string &field_name)
    {
	Teuchos::ArrayView<DataType> return_view;

	if ( field_name == "DAMPER_SOURCE_FIELD" )
	{
	    return_view = Teuchos::ArrayView<DataType>( damper->get_damping() );
	}

	return return_view;
    }

    DataType get_global_source_data(const std::string &field_name)
    {
	DataType return_val = 0.0;

	return return_val;
    }
};

} // end namespace Coupler

#endif // end damper_source_hpp
