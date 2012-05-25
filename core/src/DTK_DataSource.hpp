//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DTK_DataSource.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:43 2011
 * \brief  Interface declaration for data source applications.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATASOURCE_HPP
#define DTK_DATASOURCE_HPP

#include <string>

#include <mpi.h>

#include <Teuchos_Describable.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class DataSource
 * \brief Protocol declaration for applications acting as a data source in
 * multiphysics coupling. 
 *
 * Implementing this interface will then provide a single mesh and the fields
 * associated with that mesh. Multiple meshes will need multiple data source
 * interfaces specified.
 */
/*! 
 * \example core/test/tstInterfaces.cc
 *
 * Test of DataSource.
 */
//===========================================================================//
template<typename NodeField, typename ElementField, typename DataField>
class DataSource : Teuchos::Describable
{
  public:

    /*!
     * \brief Constructor.
     */
    DataSource()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~DataSource()
    { /* ... */ }

    /*!
     * \brief Get the communicator object for the physics implementing this
     * interface.
     * \param source_comm The MPI communicator for the source application.
     */
    virtual const MPI_Comm& getSourceComm() = 0;

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported.
     * \param field_name The name of the field for which support is being
     * checked.
     */
    virtual bool isFieldSupported( const std::string &field_name ) = 0;

    /*! 
     * \brief Provide the local source mesh nodes this includes nodes needed
     * to resolve higher order elements.
     * \param source_nodes View of the local source mesh nodes. This view
     * is required to persist. The NodeField type is expected to implement
     * FieldTraits. NodeField::value_type is expected to implement
     * NodeTraits. 
     */
    virtual const NodeField& getSourceMeshNodes() = 0;

    /*! 
     * \brief Provide the local source mesh elements.
     * \param source_elements View of the local source mesh
     * elements. This view is required to persist. The ElementField object
     * passed must be a contiguous memory view of DataTransferKit::Element
     * objects. The ElementField type is exepected to implement
     * FieldTraits. ElementField::value_type is expected to implement
     * ElementTraits. 
     */
    virtual const ElementField& getSourceMeshElements() = 0;

    /*! 
     * \brief Provide a const view of the local source data at the source mesh
     * nodes for a given field.
     * \param field_name The name of the field to provide data from.
     * \param source_node_data A persisting view of data to be sent. There are
     * two requirements for this view: 1) it is of size equal to the number of
     * source mesh nodes in the local domain, 2) the data is in the same order
     * as the nodes found provided by getSourceMeshNodes. This view is
     * required to persist. The DataField type is expected to implement
     * FieldTraits. ElementField::value_type is expected to implement
     * Teuchos::ScalarTraits.  
     */
    virtual const DataField& 
    getSourceNodeData( const std::string &field_name ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATASOURCE_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataSource.hpp
//---------------------------------------------------------------------------//
