//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Mapper.hh
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:31:05 2011
 * \brief  Mapper class defintion.
 */
//---------------------------------------------------------------------------//

#ifndef core_Mapper_hh
#define core_Mapper_hh

#include <string>
#include <vector>

#include <Mesh_Point.hpp>

#include "Transfer_Data_Field.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Tpetra_Map.hpp"

namespace coupler
{

//===========================================================================//
/*!
 * \class Mapper
 * \brief A mapper operator to generate a one way mapping for transfer between
 * a source physics and a target physics.
 *
 * \sa Mapper.cc for detailed descriptions.
 */
/*! 
 * \example core/test/tstMapper.cc
 *
 * Test of Mapper.
 */
//===========================================================================//

template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Mapper 
{
  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                               DataType;
    typedef HandleType_T                             HandleType;
    typedef CoordinateType_T                         CoordinateType;
    typedef int                                      OrdinalType;
    typedef mesh::Point<HandleType,CoordinateType>   PointType;
    typedef Transfer_Data_Field<DataType,HandleType,CoordinateType> 
                                                     Transfer_Data_Field_t;
    typedef Teuchos::RCP<Transfer_Data_Field_t>      RCP_Transfer_Data_Field;
    typedef Tpetra::Map<int>                         Tpetra_Map_t;
    typedef Teuchos::RCP<const Tpetra_Map_t>         RCP_Tpetra_Map;
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    //@}

  public:

    // Constructor.
    Mapper();

    // Destructor.
    ~Mapper();

    // Map the field from the data source onto the data target.
    void map(const RCP_Communicator comm_global,
	     RCP_Transfer_Data_Field transfer_data_field);
};

} // end namespace coupler

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "Mapper_Tpetra.t.hh"

#endif // core_Mapper_hh

//---------------------------------------------------------------------------//
//              end of core/Mapper.hh
//---------------------------------------------------------------------------//
