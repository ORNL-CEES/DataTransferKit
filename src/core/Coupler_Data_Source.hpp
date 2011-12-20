//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Coupler_Data_Source.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:43 2011
 * \brief  Interface definition for data source applications.
 */
//---------------------------------------------------------------------------//

#ifndef core_Coupler_Data_Source_hpp
#define core_Coupler_Data_Source_hpp

#include <string>

#include <Mesh_Point.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ArrayView.hpp"

namespace Coupler
{

//===========================================================================//
/*!
 * \class Data_Source
 * \brief Protocol definition for applications acting as a data source in
 * multiphysics coupling. 
 *
 * This interface is templated on the type of field data being
 * transferred, the handle type for mesh entities, and the coordinate
 * type. Oridnal type for communication is fixed to int. 
 */
/*! 
 * \example core/test/tstInterfaces.cc
 *
 * Test of Data_Source.
 */
//===========================================================================//
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Data_Source 
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
    Data_Source()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Data_Source()
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
     * \brief Given (x,y,z) coordinates and an associated globally unique
     * handle encapsulated in a point object, return true if the point is in
     * the local domain, false if not. 
     * \param point Point to query the local domain with.
     */
    virtual bool get_points(const PointType &point) = 0;

    /*! 
     * \brief Send the field data.
     * \param field_name The name of the field to send data from.
     * \return A const view of data to be sent.
     */
    virtual const Teuchos::ArrayView<DataType> 
    send_data(const std::string &field_name) = 0;

    /*!
     * \brief Given a field, set a global data element to be be sent to a
     * target.
     * \param field_name The name of the field to send data from.
     * \return The global data element.
     */
    virtual DataType set_global_data(const std::string &field_name) = 0;
};

} // end namespace Coupler

#endif // core_Coupler_Data_Source_hpp

//---------------------------------------------------------------------------//
//              end of core/Coupler_Data_Source.hpp
//---------------------------------------------------------------------------//
