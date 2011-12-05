//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Data_Field.hh
 * \author Stuart Slattery
 * \date   Fri Nov 18 11:57:58 2011
 * \brief  Transfer_Data_Field class definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Transfer_Data_Field_hh
#define core_Transfer_Data_Field_hh

#include <string>

#include "Transfer_Data_Source.hh"
#include "Transfer_Data_Target.hh"
#include "Transfer_Map.hh"

#include "Teuchos_RCP.hpp"

namespace coupler
{

//===========================================================================//
/*!
 * \class Transfer_Data_Field
 * \brief Field type for data transfers. This exists for more a explicit
 * definition of fields in the coupler user interface.
 *
 * The Transfer_Data_Field encapsulates the relationship between a
 * Transfer_Data_Source and a Transfer_Data_Target. In addition to containing
 * the map for transfer between the source and the target, the field also
 * contains indicators for the status of the mapping (transfer cannot occur
 * unless a field has been mapped) and for whether a field is distributed
 * (requiring a parallel mapping) or scalar.
 */
/*! 
 * \example core/test/tstTransfer_Data_Field.cc
 *
 * Test of Transfer_Data_Field.
 */
//===========================================================================//

template<class DataType_T>
class Transfer_Data_Field 
{

  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                               DataType;
    typedef Transfer_Data_Source<DataType>           Transfer_Data_Source_t;
    typedef Teuchos::RCP<Transfer_Data_Source_t>     RCP_Transfer_Data_Source;
    typedef Transfer_Data_Target<DataType>           Transfer_Data_Target_t;
    typedef Teuchos::RCP<Transfer_Data_Target_t>     RCP_Transfer_Data_Target;
    typedef Teuchos::RCP<Transfer_Map>               RCP_Transfer_Map;
    //@}

  private:

    // Field name.
    std::string d_field_name;

    // Data transfer source implemenation.
    RCP_Transfer_Data_Source d_source;

    // Data transfer target implemenation.
    RCP_Transfer_Data_Target d_target;

    // Topology map for transfer from the source to the target.
    RCP_Transfer_Map d_map;

    // Boolean for scalar field.
    bool d_scalar;

    // Boolean for field mapping. True if mapping complete.
    bool d_mapped;

  public:

    // Constructor.
    Transfer_Data_Field(const std::string &field_name,
			RCP_Transfer_Data_Source source,
			RCP_Transfer_Data_Target target,
			bool scalar = false);

    // Destructor.
    ~Transfer_Data_Field();

    //! Get the field name. 
    const std::string& name() 
    { return d_field_name; }

    //! Get the transfer data source.
    RCP_Transfer_Data_Source source() 
    { return d_source; }

    //! Get the transfer data target.
    RCP_Transfer_Data_Target target() 
    { return d_target; }
    
    // Set the topology map.
    void set_map(RCP_Transfer_Map transfer_map);

    //! Get the topology map.
    RCP_Transfer_Map get_map() 
    { return d_map; }

    //! Return the scalar boolean.
    bool is_scalar()
    { return d_scalar; }

    //! Return the mapped boolean.
    bool is_mapped()
    { return d_mapped; }
};

} // end namespace coupler

#endif // core_Transfer_Data_Field_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Data_Field.hh
//---------------------------------------------------------------------------//
