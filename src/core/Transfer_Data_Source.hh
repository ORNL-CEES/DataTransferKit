//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Data_Source.hh
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:43 2011
 * \brief  Interface definition for data source applications.
 */
//---------------------------------------------------------------------------//

#ifndef core_Transfer_Data_Source_hh
#define core_Transfer_Data_Source_hh

#include <vector>
#include <string>

#include <Mesh_Point.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ArrayViewDecl.hpp"

namespace coupler
{

//===========================================================================//
/*!
 * \class Transfer_Data_Source
 * \brief Definition of the interface for applications acting as a data source
 * in multiphysics coupling.
 *
 * This interface is templated on the type of field data being
 * transferred, the handle type for mesh entities, and the coordinate type.
 */
/*! 
 * \example core/test/tstTransfer_Data_Source.cc
 *
 * Test of Transfer_Data_Source.
 */
//===========================================================================//
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Transfer_Data_Source 
{
  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                                 DataType;
    typedef HandleType_T                               HandleType;
    typedef CoordinateType_T                           CoordinateType;
    typedef int                                        OrdinalType;
    typedef mesh::Point<HandleType,CoordinateType>     PointType;
    typedef Teuchos::Comm<OrdinalType>                 Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>         RCP_Communicator;
    //@}

    /*!
     * \brief Constructor.
     */
    Transfer_Data_Source()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Transfer_Data_Source()
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
     * handle, return true if the point is in the local domain, false if not.
     * \param point Point.
     */
    virtual bool get_points(PointType &point) = 0;

    /*! 
     * \brief Given an entity handle, send the field data associated with that
     * handle. 
     * \param field_name The name of the field to send data from.
     * \param points A view of the points for the data being sent.
     * \return A view of data being sent.
     */
    virtual Teuchos::ArrayView<DataType> 
    send_data(const std::string &field_name,
	      const Teuchos::ArrayView<PointType> &points) = 0;

    /*!
     * \brief Given a field, set a global data element to be be sent to a
     * target.
     * \param field_name The name of the field to send data from.
     * \return The global data element.
     */
    virtual DataType set_global_data(const std::string &field_name) = 0;
};

} // end namespace coupler

#endif // core_Transfer_Data_Source_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Data_Source.hh
//---------------------------------------------------------------------------//
