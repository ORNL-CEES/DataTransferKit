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

#ifndef coupler_Transfer_Evaluator_hh
#define coupler_Transfer_Evaluator_hh

#include "comm/global.hh"
#include <vector>
#include <string>

namespace coupler
{

//===========================================================================//
/*!
 * \class Transfer_Evaluator
 * \brief Base class definition of the transfer evaluator interface.
 *
 * This interface is templated on the type of data being transfered.
 */
//===========================================================================//
template<class FieldType_T>
class Transfer_Evaluator 
{
  public:

    //@{
    //! Useful typedefs.
    typedef FieldType_T                                   FieldType;
    typedef FieldType::value_type                         DataType;
    typedef nemesis::Communicator_t                       Communicator_t;
    typedef long int                                      Handle;
    typedef std::vector<Handle>                           Vec_Handle;
    //@}

    //! Constructor.
    Transfer_Evaluator()
    { /* ... */ }

    //! Destructor.
    virtual ~Transfer_Evaluator()
    { /* ... */ }

    //! Register communicator object.
    virtual void register_comm(Communicator_t &comm) = 0;

    //! Register a field associated with the entities. Return false if this
    //! field is not supported.
    virtual bool register_field(std::string field_name) = 0;

    //! Register coordinates with a field.
    virtual void register_xyz(std::string field_name,
			      std::vector<double> &points,
			      Vec_Handle &handles) = 0;

    //! Given (x,y,z) coordinates, return true if in the local domain, false
    //! if not. Provide a handle to the entity it was located in.
    virtual bool find_xyz(double x, 
			  double y, 
			  double z,
			  Handle &handle) = 0;

    //! Given an entity handle, get the field data associated with that
    //! handle.
    virtual void pull_data(std::string field_name,
			   Handle handle,
			   DataType &data);

    //! Given an entity handle, set the field data associated with that
    //! handle.
    virtual void push_data(std::string field_name,
			   Handle handle, 
			   DataType data);

    //! Perfom a global integration on a field for rebalance.
    virtual void integrate_field(std::string field_name,
				 DataType &field_norm) = 0;

    //! Perform a rebalance on a field for conservation.
    virtual void rebalance(std::string field_name,
			   DataType field_norm) = 0;
};

} // end namespace coupler

#endif // coupler_Transfer_Evaluator_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Evaluator.hh
//---------------------------------------------------------------------------//
