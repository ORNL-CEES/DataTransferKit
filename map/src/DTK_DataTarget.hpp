//---------------------------------------------------------------------------//
/*!
 * \file   DTK_DataTarget.hpp
 * \author Stuart R. Slattery
 * \brief  DataTarget declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATATARGET_HPP
#define DTK_DATATARGET_HPP

#include <string>
#include <map>

#include <mpi.h>

#include "DTK_FieldTraits.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class DataTarget
 * \brief DataTarget objects provide coordinate fields on which to evaluate
 * data fields and memory spaces to write the evaluated data.
 */
//===========================================================================//
template<class CoordinateField, class DataField>
class DataTarget
{
  public:

    //@{
    //! Typedefs.
    typedef DataField                        data_field_type;
    typedef FieldTraits<DataField>           DataFT;
    typedef CoordinateField                  coordinate_field_type;
    typedef FieldTraits<CoordinateField>     CoordinateFT;
    typedef Teuchos::Comm<int>               CommType;
    typedef Teuchos::RCP<CommType>           RCP_Comm;
    //@}

  public:

    // Constructor.
    DataTarget( const CoordinateField& coordinate_field,
		const MPI_Comm& mpi_comm );
		
    // Destructor.
    ~DataTarget();

    // Register a target field.
    void registerTargetField( const std::string& field_name, 
			      DataField& data_space );

    //! Get the id for a field given its name.
    std::size_t getFieldId( const std::string& name ) const
    { return d_name_map.find( name )->second; }

    //! Get the data space for a field given its id.
    DataField& getDataSpace( const std::size_t id ) const
    { return d_data_map.find( id )->second; }

    //! Get the coordinate field.
    const CoordinateField& getCoordinateField() const
    { return d_coord_field; }

    //! Get the communicator.
    const RCP_Comm& getComm() const
    { return d_comm; }

  private:

    // Coordinate field.
    CoordinateField d_coord_field;

    // Communicator.
    RCP_Comm d_comm;

    // Name to id map.
    std::map<std::string,std::size_t> d_name_map;

    // Id to data field map.
    std::map<std::size_t,DataField> d_data_map;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_DataTarget_def.hpp"

#endif // DTK_DATATARGET_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataTarget.hpp
//---------------------------------------------------------------------------//
