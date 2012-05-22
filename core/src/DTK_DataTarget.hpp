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

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>
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
template<class DataType, class HandleType, class CoordinateType, int DIM>
class DataTarget : public Teuchos::Describable
{
  public:

    //@{
    //! Useful typedefs.
    typedef int                                        OrdinalType;
    typedef Node<DIM,HandleType,CoordinateType>        NodeType;
    typedef Teuchos::Comm<OrdinalType>                 Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>         RCP_Communicator;
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
     * \return The communicator for this physics.
     */
    virtual RCP_Communicator getTargetComm() = 0;

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
     * returned from the transfer operation. This view required to persist.
     * registered with.
     * \return View of the local target nodes.
     */
    virtual const Teuchos::ArrayRCP<NodeType> 
    getTargetMeshNodes() = 0;

    /*! 
     * \brief Provide a persisting, non-const view of the local data vector 
     * associated with the nodes provided by get_target_nodes.
     * \param field_name The name of the field to receive data from.
     * \return A non-const persisting view of the data vector to be
     * populated. This view has two requirements: 1) It is of size equal to
     * the number of nodes provided by get_target_nodes, 2) It is a
     * persisting view that will be used to write data into the underlying
     * vector. The order of the data provided will be in the same order as the
     * local nodes provided by get_target_nodes.
     */
    virtual Teuchos::ArrayRCP<DataType> 
    getTargetDataSpace( const std::string &field_name ) = 0;

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
