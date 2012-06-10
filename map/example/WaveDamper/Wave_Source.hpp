#ifndef wave_source_hpp
#define wave_source_hpp

#include <algorithm>

#include "Wave.hpp"

#include <DataTransferKit_Point.hpp>
#include <DataTransferKit_DataSource.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"

//---------------------------------------------------------------------------//
namespace DataTransferKit {

// DataSource interface implementation for the Wave code.
template<class DataType, class HandleType, class CoordinateType, int DIM>
class Wave_DataSource 
    : public DataSource<DataType,HandleType,CoordinateType,DIM>
{
  public:

    typedef int                                      OrdinalType;
    typedef Point<1,HandleType,CoordinateType>         PointType;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    typedef Teuchos::RCP<Wave>                       RCP_Wave;

  private:

    // Wave object to operate on.
    RCP_Wave wave;

  public:

    Wave_DataSource(RCP_Wave _wave)
	: wave(_wave)
    { /* ... */ }

    ~Wave_DataSource()
    { /* ... */ }

    RCP_Communicator get_source_comm()
    {
	return wave->get_comm();
    }

    bool is_field_supported(const std::string &field_name)
    {
	bool return_val = false;

	if (field_name == "WAVE_SOURCE_FIELD")
	{
	    return_val = true;
	}

	return return_val;
    }

    bool is_local_point(const PointType &test_point)
    {
	bool return_val = false;

	if ( std::find(wave->get_grid().begin(), 
		       wave->get_grid().end(),
		       test_point.getCoords()[0] )
	     != wave->get_grid().end() )
	{
	    return_val = true;
	}

	return return_val;
    }

    const Teuchos::ArrayRCP<DataType> 
    get_source_data(const std::string &field_name)
    {
	Teuchos::ArrayRCP<DataType> return_view;

	if ( field_name == "WAVE_SOURCE_FIELD" )
	{
	    return_view = Teuchos::arcp<DataType>( 
		&wave->get_f()[0], 0, (int) wave->get_f().size(), false  );
	}

	return return_view;
    }

    DataType get_global_source_data(const std::string &field_name)
    {
	DataType return_val = 0.0;

	return return_val;
    }
};

} // end namespace DataTransferKit

#endif // end wave_source_hpp
