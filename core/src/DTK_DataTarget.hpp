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

#include "DTK_FieldTraits.hpp"

#include <Teuchos_Describable.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class DataTarget
 * \brief Protocol definition for applications acting as a data target in
 * multiphysics coupling. DTK will neither allocate or deallocate memory in
 * the target application throught this interface.
 *
 */
/*! 
 * \example test/tstInterfaces.cc
 *
 * Test of DataTarget.
 */
//===========================================================================//
template<typename CoordinateField, typename DataField>
class DataTarget : public Teuchos::Describable
{
  public:

    //@{
    //! Coordinate field typedefs.
    typedef CoordinateField                              coordinate_field_type;
    typedef typename FieldTraits<CoordinateField>::value_type  coordinate_type;
    //@}

    //@{
    //! Data typedefs.
    typedef DataField                                      data_field_type;
    typedef typename FieldTraits<DataField>::value_type    data_type;
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
     * \brief Check whether or not a field is supported.
     * \param field_name The name of the field for which support is being
     * checked.
     * \return Return false if this field is not supported. 
     */
    virtual bool isFieldSupported( const std::string &field_name ) = 0;

    /*!
     * \brief Provide the target mesh node coordinates to which data will be
     * transferred.
     * The order of these coordinates will correspond to the order of the data
     * returned from the transfer operation.  
     * \param target_nodes View of the local target nodes. This view required
     * to persist. The CoordinateField type is expected to implement
     * FieldTraits. CoordinateField::value_type is required to be of type
     * double. The coordinates are required to be three dimensional. The
     * coordinates are expected to be interleaved. 
     * ( x0, y0, z0,x1, y1, z1, ... , xN, yN, zN )
     */
    virtual const CoordinateField& getTargetCoordinates() = 0;

    /*! 
     * \brief Provide a persisting, non-const view of the local data vector
     * associated with the nodes provided by getTargetCoordinates().
     * \param field_name The name of the field to receive data from. 
     * \return A non-const persisting view of the data vector to be
     * populated. This view has two requirements: 1) It is of size equal to
     * the number of nodes provided by getTargetCoordinates(), 2) It is a
     * persisting view that will be used to write data into the underlying
     * vector. The order of the data provided will be in the same order as the
     * local coordinates provided by getTargetCoordinates(). The DataField type is
     * expected to implement FieldTraits. DataField::value_type is expected to
     * implement Teuchos::ScalarTraits.
     */
    virtual DataField& getTargetDataSpace( const std::string &field_name ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATATARGET_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataTarget.hpp
//---------------------------------------------------------------------------//
