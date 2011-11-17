//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Data_Source.hh
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:43 2011
 * \brief  Interface definition for data source applications.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Transfer_Data_Source_hh
#define core_Transfer_Data_Source_hh

#include <vector>
#include <string>

#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Transfer_Data_Source
 * \brief Definition of the interface for applications acting as a data source
 * in multiphysics coupling.
 *
 * This interface is templated on the type of field data being
 * transfered. Handle type is fixed to integer and coordinate type is fixed to
 * double. These could be templated in the future.
 */
/*! 
 * \example core/test/tstTransfer_Data_Source.cc
 *
 * Test of Transfer_Data_Source.
 */
//===========================================================================//
template<class DataType_T>
class Transfer_Data_Source 
{
  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                               DataType;
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      HandleType;
    typedef double                                   CoordinateType;
    //@}

    /*!
     * \brief Constructor.
     */
    Data_Transfer_Source()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Data_Transfer_Source()
    { /* ... */ }

    /*!
     * \brief Register communicator object.
     * \param comm The communicator for this physics.
     */
    virtual void register_comm(const Communicator &comm) = 0;

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
     * \param handle The globally unique handle associated with the point.
     * \param x X coordinate.
     * \param y Y coordinate.
     * \param z Z coordinate.
     * coordinates.
     */
    virtual bool get_points(HandleType handle,
			    CoordinateType x, 
			    CoordinateType y,
			    CoordinateType z) = 0;

    /*! 
     * \brief Given an entity handle, send the field data associated with that
     * handle. 
     * \param field_name The name of the field to send data from.
     * \param handles The enitity handles for the data being sent.
     * \param data The data being sent.
     */
    virtual void send_data(const std::string &field_name,
			   const std::vector<HandleType> &handles,
			   std::vector<DataType> &data) = 0;

    /*!
     * \brief Given a field, set a global data element to be be sent to a
     * target.
     * \param field_name The name of the field to send data from.
     * \param data The global data element.
     */
    virtual void set_global_data(const std::string &field_name,
				 DataType &data) = 0;
};

} // end namespace coupler

#endif // core_Transfer_Data_Source_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Data_Source.hh
//---------------------------------------------------------------------------//
