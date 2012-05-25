//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DTK_DataTarget.hpp
 * \author Stuart R. Slattery
 * \date   Thu Nov 17 07:53:54 2011
 * \brief  Interface definition for transfer data targets.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATATARGET_HPP
#define DTK_DATATARGET_HPP

#include <string>

#include <mpi.h>

#include "FieldTraits.hpp"

#include <Teuchos_Describable.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class DataTarget
 * \brief Protocol definition for applications acting as a data target in
 * multiphysics coupling. 
 *
 */
/*! 
 * \example test/tstInterfaces.cc
 *
 * Test of DataTarget.
 */
//===========================================================================//
template<typename NodeField, typename DataField>
class DataTarget : public Teuchos::Describable
{
  public:

    //@{
    //! Typedefs.
    typedef NodeField                             node_field_type;
    typedef FieldTraits<NodeField>::value_type    node_type;
    typedef DataField                             data_field_type;
    typedef FieldTraits<DataField>::value_type    data_type;
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
     * \param target_comm The MPI communicator for the target application.
     */
    virtual const MPI_Comm& getTargetComm() = 0;

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
     * returned from the transfer operation.  \param target_nodes View of the
     * local target nodes. This view required to persist. The NodeField
     * type is expected to implement FieldTraits. NodeField::value_type
     * is expected to implement NodeTraits. 
     */
    virtual const NodeField& getTargetNodes() = 0;

    /*! 
     * \brief Provide a persisting, non-const view of the local data vector
     * associated with the nodes provided by get_target_nodes.  
     * \param field_name The name of the field to receive data from. 
     * \return A non-const persisting view of the data vector to be
     * populated. This view has two requirements: 1) It is of size equal to
     * the number of nodes provided by getTargetNodes(), 2) It is a persisting
     * view that will be used to write data into the underlying vector. The
     * order of the data provided will be in the same order as the local nodes
     * provided by getTargetNodes(). The DataField type is expected to
     * implement FieldTraits. ElementField::value_type is expected to
     * implement Teuchos::ScalarTraits.
     */
    virtual DataField& getTargetDataSpace( const std::string &field_name ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATATARGET_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataTarget.hpp
//---------------------------------------------------------------------------//
