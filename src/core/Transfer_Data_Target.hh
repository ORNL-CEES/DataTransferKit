//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Data_Target.hh
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:54 2011
 * \brief  Interface definition for transfer data targets.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Transfer_Data_Target_hh
#define core_Transfer_Data_Target_hh

namespace coupler
{

//===========================================================================//
/*!
 * \class Transfer_Data_Target
 * \brief Interface defintion for transfer data targets.
 *
 * This interface is templated on the type of field data being
 * transfered. Handle type is fixed to integer and coordinate type is fixed to
 * double. These could be templated in the future.
 */
/*! 
 * \example core/test/tstTransfer_Data_Target.cc
 *
 * Test of Transfer_Data_Target.
 */
//===========================================================================//
template<class DataType_T>
class Transfer_Data_Target 
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
    Data_Transfer_Target()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Data_Transfer_Target()
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
     * \brief Set cartesian coordinates with a field. The coordinate
     * vector should be interleaved. The handle vector should consist of
     * globally unique handles. 
     * \param field_name The name of the field that the coordinates are being
     * registered with.
     * \param handles Point handle array.
     * \param coordinates Point coordinate array.
     */
    virtual void set_points(const std::string &field_name,
			    std::vector<HandleType> &handles,
			    std::vector<CoordinateType> &coordinates) = 0;

    /*! 
     * \brief Given an entity handle, send the field data associated with that
     * handle. 
     * \param field_name The name of the field to send data from.
     * \param handles The enitity handles for the data being sent.
     * \param data The data being sent.
     */
    virtual void receive_data(const std::string &field_name,
			      const std::vector<HandleType> &handles,
			      const std::vector<DataType> &data) = 0;

    /*!
     * \brief Given a field, get a global data element to be be sent to a
     * target.
     * \param field_name The name of the field to send data from.
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
