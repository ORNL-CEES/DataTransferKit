//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DataTransferKit_DataTarget.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:54 2011
 * \brief  Interface definition for transfer data targets.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATATARGET_HPP
#define DTK_DATATARGET_HPP

#include <string>

#include "DataTransferKit_Point.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
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
    typedef Point<DIM,HandleType,CoordinateType>       PointType;
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
    virtual RCP_Communicator get_target_comm() = 0;

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported. 
     * \param field_name The name of the field for which support is being
     * checked.
     */
    virtual bool is_field_supported( const std::string &field_name ) = 0;

    /*!
     * \brief Given a field, provide the local points to map data on to. The
     * order of these points will correspond to the order of the data returned
     * from the transfer operation. This view required to persist.
     * \param field_name The name of the field that the points are being
     * registered with.
     * \return View of the local target points.
     */
    virtual const Teuchos::ArrayRCP<PointType> 
    get_target_points( const std::string &field_name ) = 0;

    /*! 
     * \brief Provide a persisting, non-const view of the local data vector 
     * associated with the points provided by get_target_points.
     * \param field_name The name of the field to receive data from.
     * \return A non-const persisting view of the data vector to be
     * populated. This view has two requirements: 1) It is of size equal to
     * the number of points provided by get_target_points, 2) It is a
     * persisting view that will be used to write data into the underlying
     * vector. The order of the data provided will be in the same order as the
     * local points provided by get_target_points. It is reccomended that the
     * return ArrayRCP not own the underlying memory as the destructors for
     * other objects used by DataTransferKit (i.e. Tpetra::vector) will
     * attempt to deallocate it.
     */
    virtual Teuchos::ArrayRCP<DataType> 
    get_target_data_space( const std::string &field_name ) = 0;

    /*!
     * \brief Given a field, set a global data element provided by the
     * source.
     * \param field_name The name of the field to set data to.
     * \param data The provided global data element.
     */
    virtual void set_global_target_data( const std::string &field_name,
					 const DataType &data ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATATARGET_HPP

//---------------------------------------------------------------------------//
//              end of DataTransferKit_DataTarget.hpp
//---------------------------------------------------------------------------//
