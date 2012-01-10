//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Coupler_DataTarget.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:54 2011
 * \brief  Interface definition for transfer data targets.
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_DATATARGET_HPP
#define COUPLER_DATATARGET_HPP

#include <string>

#include "Coupler_Point.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Coupler
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
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class DataTarget 
{
  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                                 DataType;
    typedef HandleType_T                               HandleType;
    typedef CoordinateType_T                           CoordinateType;
    typedef int                                        OrdinalType;
    typedef Point<HandleType,CoordinateType>           PointType;
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
    virtual bool is_field_supported(const std::string &field_name) = 0;

    /*!
     * \brief Given a field, provide the local points to map data on to. The
     * order of these points will correspond to the order of the data returned
     * from the transfer operation. This view is not required to persist.
     * \param field_name The name of the field that the points are being
     * registered with.
     * \return View of the local target points.
     */
    virtual const Teuchos::ArrayView<PointType> 
    get_target_points(const std::string &field_name) = 0;

    /*! 
     * \brief Provide a persisting, non-const view of the local data vector 
     * associated with the points provided by get_target_points.
     * \param field_name The name of the field to receive data from.
     * \return A non-const view of the data vector to be populated. This view
     * has two requirements: 1) It is of size equal to the number of points
     * provided by get_target_points, 2) It is a persistent view that will be
     * used to write data into the underlying vector. The order of the data
     * provided will be in the same order as the local points provided by
     * get_target_points. 
     */
    virtual Teuchos::ArrayView<DataType> 
    get_target_data_space(const std::string &field_name) = 0;

    /*!
     * \brief Given a field, set a global data element provided by the
     * source.
     * \param field_name The name of the field to set data to.
     * \param data The provided global data element.
     */
    virtual void set_global_target_data(const std::string &field_name,
					const DataType &data) = 0;
};

} // end namespace Coupler

#endif // COUPLER_DATATARGET_HPP

//---------------------------------------------------------------------------//
//              end of Coupler_DataTarget.hpp
//---------------------------------------------------------------------------//
