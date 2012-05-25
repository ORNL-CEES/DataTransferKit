//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DTK_DataSource.hpp
 * \author Stuart R. Slattery
 * \date   Thu Nov 17 07:53:43 2011
 * \brief  Interface declaration for data source applications.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATASOURCE_HPP
#define DTK_DATASOURCE_HPP

#include <string>
#include <vector>

#include <mpi.h>

#include "DTK_NodeTraits.hpp"
#include "DTK_ElementTraits.hpp"
#include "DTK_FieldTraits.hpp"

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

    //@{
    //! Typdefs.
    typedef NodeField                                         node_field_type;
    typedef typename FieldTraits<NodeField>::value_type       node_type;
    typedef typename NodeTraits<node_type>::coordinate_type   node_coordinate_type;
    typedef ElementField                                      element_field_type;
    typedef typename FieldTraits<ElementField>::value_type    element_type;
    typedef typename ElementTraits<element_type>::handle_type element_handle_type;
    typedef DataField                                         data_field_type;
    typedef typename FieldTraits<DataField>::value_type       data_type;
    typedef std::vector<node_coordinate_type>                 CoordinateVector;
    typedef std::vector<element_handle_type>                  HandleVector;
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
     * is required to persist. The NodeField type is expected to
     * implement FieldTraits. NodeField::value_type is expected to
     * implement NodeTraits. 
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
     * \brief Provide a copy of the local source data evaluated in the
     * given source mesh elements at the given nodes.
     * \param field_name The name of the field to provide data from.
     * \param element_handles The local element handles in which to evaluate
     * the field. See typedef for type construction.
     * \param node_field The nodes at which to evaluate the field. This vector
     * will be the length of the element handle vector times the
     * dimensionality of the coordinate. See typedef for type construction.
     * \return The evaluated field data at the nodes. The DataField type
     * is expected to implement FieldTraits. DataField::value_type is
     * expected to implement Teuchos::ScalarTraits. 
     */
    virtual const DataField
    evaluateFieldOnTargetNodes( const std::string& field_name,
				const HandleVector& element_handles,
				const CoordinateVector& node_coordinates ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATASOURCE_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataSource.hpp
//---------------------------------------------------------------------------//
