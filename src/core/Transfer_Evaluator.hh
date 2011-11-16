//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Evaluator.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:59 2011
 * \brief  Transfer_Evaluator class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Transfer_Evaluator_hh
#define core_Transfer_Evaluator_hh

#include <vector>
#include <string>

#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Transfer_Evaluator
 * \brief Definition of the transfer evaluator interface for multiphysics
 * coupling.
 *
 * This interface is templated on the type of field data being
 * transfered. Handle type is fixed to integer and coordinate type is fixed to
 * double. These could be templated in the future.
 */
//===========================================================================//
template<class DataType_T>
class Transfer_Evaluator 
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
    Transfer_Evaluator()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Transfer_Evaluator()
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
     * \brief Register cartesian coordinates with a field. The coordinate
     * vector should be interleaved. The handle vector should consist of
     * globally unique handles. 
     * \param field_name The name of the field that the coordinates are being
     * registered with.
     * \param handles Point handle array.
     * \param coordinates Point coordinate array.
     */
    virtual void register_points(const std::string &field_name,
				 std::vector<HandleType> &handles,
				 std::vector<CoordinateType> &coordinates) = 0;

    /*! 
     * \brief Given (x,y,z) coordinates and an associated globally unique
     * handle, return true if the point is in the local domain, false if not.
     * \param handle The globally unique handle associated with the point.
     * \param x X coordinate.
     * \param y Y coordinate.
     * \param z Z coordinate.
     * coordinates.
     */
    virtual bool find_point(HandleType handle,
			    CoordinateType x, 
			    CoordinateType y,
			    CoordinateType z) = 0;

    /*! 
     * \brief Given an entity handle, get the field data associated with that
     * handle. 
     * \param field_name The name of the field to pull data from.
     * \param handle The enitity handles for the data being pulled.
     * \param data The data being pulled.
     */
    virtual void pull_data(const std::string &field_name,
			   std::vector<HandleType> &handle,
			   std::vector<DataType> &data) = 0;

    /*!
     * \brief Given an entity handle, set the field data associated with that
     * handle. 
     * \param field_name The name of the field to push data from.
     * \param handle The enitity handles for the data being pushed.
     * \param data The data being pushed.
     */
    virtual void push_data(const std::string &field_name,
			   const std::vector<HandleType> &handle, 
			   const std::vector<DataType> &data) = 0;

    /*!
     * \brief Perfom a global integration on a field for rebalance.
     * \param field_name The name of the field being balanced.
     * \param field_norm The normalization factor computed by the global
     * integration.
     */
    virtual void integrate(const std::string &field_name,
			   DataType &field_norm) = 0;

    /*!
     * \brief Perfom a global rebalance on a field given a normalization
     * factor.
     * \param field_name The name of the field being balanced.
     * \param field_norm The normalization factor computed by the global
     * integration.
     */
    virtual void rebalance(const std::string &field_name,
			   DataType field_norm) = 0;
};

} // end namespace coupler

#endif // core_Transfer_Evaluator_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Evaluator.hh
//---------------------------------------------------------------------------//
