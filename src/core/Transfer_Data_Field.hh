//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Data_Field.hh
 * \author Stuart Slattery
 * \date   Fri Nov 18 11:57:58 2011
 * \brief  Transfer_Data_Field class definition.
 */
//---------------------------------------------------------------------------//

#ifndef core_Transfer_Data_Field_hh
#define core_Transfer_Data_Field_hh

#include <string>

#include "Transfer_Data_Source.hh"
#include "Transfer_Data_Target.hh"

#include "Teuchos_RCP.hpp"

#include "Tpetra_Map.hpp"

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

template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Transfer_Data_Field
{

  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                               DataType;
    typedef HandleType_T                             HandleType;
    typedef CoordinateType_T                         CoordinateType;
    typedef Transfer_Data_Source<DataType,HandleType,CoordinateType> Transfer_Data_Source_t;
    typedef Teuchos::RCP<Transfer_Data_Source_t>     RCP_Transfer_Data_Source;
    typedef Transfer_Data_Target<DataType,HandleType,CoordinateType> Transfer_Data_Target_t;
    typedef Teuchos::RCP<Transfer_Data_Target_t>     RCP_Transfer_Data_Target;
    typedef Tpetra::Map<int>                         Tpetra_Map_t;
    typedef Teuchos::RCP<Tpetra_Map_t>               RCP_Tpetra_Map;
    //@}

  private:

    // Field name.
    std::string d_field_name;

    // Data transfer source implemenation.
    RCP_Transfer_Data_Source d_source;

    // Data transfer target implemenation.
    RCP_Transfer_Data_Target d_target;

    // Tpetra map for the data source.
    RCP_Tpetra_Map d_source_map;

    // Tpetra map for the data target.
    RCP_Tpetra_Map d_target_map;

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
    
    //! Set the maps for the data source and target.
    void set_maps(RCP_Tpetra_Map source_map, RCP_Tpetra_Map target_map);

    //! Get the map for the data source.
    RCP_Tpetra_Map get_source_map()
    { return d_source_map; }

    //! Get the map for the data target.
    RCP_Tpetra_Map get_target_map()
    { return d_target_map; }

    //! Return the scalar boolean.
    bool is_scalar()
    { return d_scalar; }

    //! Return the mapped boolean.
    bool is_mapped()
    { return d_mapped; }
};

} // end namespace coupler

//---------------------------------------------------------------------------//
// TEMPLATE MEMBERS
//---------------------------------------------------------------------------//

#include "Transfer_Data_Field.t.hh"

#endif // core_Transfer_Data_Field_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Data_Field.hh
//---------------------------------------------------------------------------//
