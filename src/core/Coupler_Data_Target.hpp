//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Coupler_Data_Target.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:54 2011
 * \brief  Interface definition for transfer data targets.
 */
//---------------------------------------------------------------------------//

#ifndef core_Coupler_Data_Target_hpp
#define core_Coupler_Data_Target_hpp

#include <string>

#include <Mesh_Point.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ArrayView.hpp"

namespace Coupler
{

//===========================================================================//
/*!
 * \class Data_Target
 * \brief Protocol definition for applications acting as a data target in
 * multiphysics coupling. 
 *
 * This interface is templated on the type of field data being
 * transferred, the handle type for mesh entities, and the coordinate
 * type. Ordinal type for communication is int.
 */
/*! 
 * \example core/test/tstInterfaces.cc
 *
 * Test of Data_Target.
 */
//===========================================================================//
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Data_Target 
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
    Data_Target()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Data_Target()
    { /* ... */ }

    /*!
     * \brief Register communicator object.
     * \return The communicator for this physics.
     */
    virtual RCP_Communicator comm() = 0;

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported. 
     * \param field_name The name of the field for which support is being
     * checked.
     */
    virtual bool field_supported(const std::string &field_name) = 0;

    /*!
     * \brief Set cartesian coordinates with a field. The order of these
     * points will correspond to the order of the data returned from the
     * transfer operation.
     * \param field_name The name of the field that the points are being
     * registered with.
     * \return View of the local target points.
     */
    virtual const Teuchos::ArrayView<PointType> 
    set_points(const std::string &field_name) = 0;

    /*! 
     * \brief Receive the field data by providing a view to be populated. 
     * \param field_name The name of the field to receive data from.
     * \return A non-const view of the data vector to be populated.
     */
    virtual Teuchos::ArrayView<DataType> 
    receive_data(const std::string &field_name) = 0;

    /*!
     * \brief Given a field, get a global data element to be be received from
     * a source.
     * \param field_name The name of the field to receive data from.
     * \param data The global data element.
     */
    virtual void get_global_data(const std::string &field_name,
				 const DataType &data) = 0;
};

} // end namespace Coupler

#endif // core_Coupler_Data_Target_hpp

//---------------------------------------------------------------------------//
//              end of core/Coupler_Data_Target.hpp
//---------------------------------------------------------------------------//
