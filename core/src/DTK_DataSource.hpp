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
template<typename SourceNodeField, typename TargetNodeField,
	 typename SourceElementField, typename SourceDataField>
class DataSource : Teuchos::Describable
{
  public:

    //@{
    //! Typdefs.
    typedef SourceNodeField              source_node_field_type;
    typedef TargetNodeField              target_node_field_type;
    typedef SourceElementField           source_element_field_type;
    typedef SourceDataField              source_data_field_type;
    //@}

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
     * is required to persist. The SourceNodeField type is expected to
     * implement FieldTraits. SourceNodeField::value_type is expected to
     * implement NodeTraits. 
     */
    virtual const SourceNodeField& getSourceMeshNodes() = 0;

    /*! 
     * \brief Provide the local source mesh elements.
     * \param source_elements View of the local source mesh
     * elements. This view is required to persist. The ElementField object
     * passed must be a contiguous memory view of DataTransferKit::Element
     * objects. The SourceElementField type is exepected to implement
     * FieldTraits. SourceElementField::value_type is expected to implement
     * ElementTraits.  
     */
    virtual const SourceElementField& getSourceMeshElements() = 0;

    /*! 
     * \brief Provide a const view of the local source data evaluated in the
     * given source mesh elements at the given nodes.
     * \param field_name The name of the field to provide data from.
     * \param element_field The local elements in which to evaluate the field.
     * \param node_field The nodes at which to evaluate the field. This type
     * can be expected to implement FieldTraits and
     * TargetNodeField::value_type can be expected to implement NodeTraits.
     * \return The evaluated field data at the nodes. The DataField type is
     * expected to implement FieldTraits. ElementField::value_type is expected
     * to implement Teuchos::ScalarTraits.  
     */
    virtual const SourceDataField&
    evaluateFieldOnTargetNodes( const std::string& field_name,
				const SourceElementField& element_field,
				const TargetNodeField& node_field ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATASOURCE_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataSource.hpp
//---------------------------------------------------------------------------//
