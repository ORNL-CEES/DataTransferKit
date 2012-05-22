//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DTK_DataTarget.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:54 2011
 * \brief  Interface definition for transfer data targets.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATATARGET_HPP
#define DTK_DATATARGET_HPP

#include <string>

#include "DTK_Node.hpp"

#include <Teuchos_Describable.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class DataTarget
 * \brief Protocol definition for applications acting as a data target in
 * multiphysics coupling. 
 *
 * This interface is templated on the type of field data being
 * transferred, the handle type for mesh entities, and the coordinate
 * type. Ordinal type for communication is int.
 */
/*! 
 * \example test/tstInterfaces.cc
 *
 * Test of DataTarget.
 */
//===========================================================================//
template<int DIM, class DataType, class HandleType, class CoordinateType>
class DataTarget : public Teuchos::Describable
{
  public:

    //@{
    //! Useful typedefs.
    typedef int                                        OrdinalType;
    typedef Node<DIM,HandleType,CoordinateType>        NodeType;
    //@}

    /*!
     * \brief Constructor.
     */
    DataTarget()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~DataTarget()
    { /* ... */ }

    /*!
     * \brief Get the communicator object for the physics implementing this
     * interface.
     * \param target_comm The communicator for the target application. This
     * class is requried to implement MPI primitives.
     */
    template<class Communicator>
    virtual void getTargetComm( const Communicator &target_comm ) = 0;

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported. 
     * \param field_name The name of the field for which support is being
     * checked.
     */
    virtual bool isFieldSupported( const std::string &field_name ) = 0;

    /*!
     * \brief Provide the target mesh nodes to which data will be transferred.
     * The order of these nodes will correspond to the order of the data
     * returned from the transfer operation. 
     * \param target_nodes View of the local target nodes. This view
     * required to persist. The NodeField object passed must be a contiguous
     * memory view of DataTransferKit::Node objects. The NodeField object must
     * have a .size() function that returns the number of nodes in the
     * contiguous memory space. The NodeField object must have a [] operator
     * for element access.
     */
    template<class NodeField>
    virtual void getTargetMeshNodes( const NodeField& target_nodes ) = 0;

    /*! 
     * \brief Provide a persisting, non-const view of the local data vector 
     * associated with the nodes provided by get_target_nodes.
     * \param field_name The name of the field to receive data from.
     * \return A non-const persisting view of the data vector to be
     * populated. This view has two requirements: 1) It is of size equal to
     * the number of nodes provided by get_target_nodes, 2) It is a
     * persisting view that will be used to write data into the underlying
     * vector. The order of the data provided will be in the same order as the
     * local nodes provided by get_target_nodes. The DataField object passed
     * must be a contiguous memory view of DataTransferKit::Data objects. The
     * DataField object must have a .size() function that returns the number
     * of datas in the contiguous memory space. The DataField object must have
     * a [] operator for element access.
     */
    template<class DataField>
    virtual void getTargetDataSpace( const std::string &field_name,
				     DataField &target_data_space ) = 0;

    /*!
     * \brief Given a field, set a global data element provided by the
     * source.
     * \param field_name The name of the field to set data to.
     * \param data The provided global data element.
     */
    virtual void setGlobalTargetData( const std::string &field_name,
				      const DataType &data ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATATARGET_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataTarget.hpp
//---------------------------------------------------------------------------//
