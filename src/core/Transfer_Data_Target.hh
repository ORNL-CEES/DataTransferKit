//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Data_Target.hh
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:54 2011
 * \brief  Interface definition for transfer data targets.
 */
//---------------------------------------------------------------------------//

#ifndef core_Transfer_Data_Target_hh
#define core_Transfer_Data_Target_hh

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
 * \class Transfer_Data_Target
 * \brief Definition of the interface for applications acting as a data target
 * in multiphysics coupling.
 *
 * This interface is templated on the type of field data being
 * transferred, the handle type for mesh entities, and the coordinate type.
 */
/*! 
 * \example core/test/tstInterfaces.cc
 *
 * Test of Transfer_Data_Target.
 */
//===========================================================================//
template<class DataType_T, class HandleType_T, class CoordinateType_T>
class Transfer_Data_Target 
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
    Transfer_Data_Target()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Transfer_Data_Target()
    { /* ... */ }

    /*!
     * \brief Register communicator object.
     * \return The communicator for this physics.
     */
    virtual const RCP_Communicator comm() = 0;

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported. 
     * \param field_name The name of the field for which support is being
     * checked.
     */
    virtual bool field_supported(const std::string &field_name) = 0;

    /*!
     * \brief Set cartesian coordinates with a field. The coordinate
     * vector should be interleaved. The handle vector should consist of
     * globally unique handles. 
     * \param field_name The name of the field that the coordinates are being
     * registered with.
     * \return Array view into a point array.
     */
    virtual const Teuchos::ArrayView<PointType> 
    set_points(const std::string &field_name) = 0;

    /*! 
     * \brief Given an entity handle, receive the field data associated with
     * that handle. 
     * \param field_name The name of the field to receive data from.
     * \param handles The enitity handles for the data being received.
     * \param data The data being received.
     */
    virtual void receive_data(const std::string &field_name,
			      const Teuchos::ArrayView<PointType> &points,
			      const Teuchos::ArrayView<DataType> &data) = 0;

    /*!
     * \brief Given a field, get a global data element to be be received from
     * a source.
     * \param field_name The name of the field to receive data from.
     * \param data The global data element.
     */
    virtual void get_global_data(const std::string &field_name,
				 const DataType &data) = 0;
};

} // end namespace coupler

#endif // core_Transfer_Data_Target_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Data_Target.hh
//---------------------------------------------------------------------------//
