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
 * This interface is templated on the type of field data being
 * transfered. Handle type is fixed to integer and coordinate type is fixed to
 * double. These could be templated in the future.
 */
//===========================================================================//
template<class FieldType_T>
class Transfer_Evaluator 
{
  public:

    //@{
    //! Useful typedefs.
    typedef FieldType_T                              FieldType;
    typedef nemesis::Communicator_t                  Communicator;
    typedef int                                      Handle;
    typedef const Handle*                            Handle_Iterator;
    typedef double                                   Coordinate;
    typedef const Coordinate*                        Coord_Iterator;
    //@}

    //! Constructor.
    Transfer_Evaluator()
    { /* ... */ }

    //! Destructor.
    virtual ~Transfer_Evaluator()
    { /* ... */ }

    //! Register communicator object.
    virtual void register_comm(Communicator &comm) = 0;

    //! Register a field associated with the entities. Return false if this
    //! field is not supported.
    virtual bool register_field(std::string field_name) = 0;

    //! Register cartesian coordinates with a field. The coordinate vector
    //! should be interleaved. The handle vector should consist of globally
    //! unique handles. These iterators imply contiguous memory storage.
    virtual void register_xyz(std::string field_name,
			      Coord_Iterator &points_begin,
			      Coord_Iterator &points_end,
			      Handle_Iterator &handles_begin,
			      Handle_Iterator &handles_end) = 0;

    //! Given (x,y,z) coordinates and an associated globally unique handle,
    //! return true if in the local domain, false if not.
    virtual bool find_xyz(Coordinate x, 
			  Coordinate y,
			  Coordinate z,
			  Handle handle) = 0;

    //! Given an entity handle, get the field data associated with that
    //! handle.
    virtual void pull_data(std::string field_name,
			   Handle handle,
			   FieldType &data) = 0;

    //! Given an entity handle, set the field data associated with that
    //! handle.
    virtual void push_data(std::string field_name,
			   Handle handle, 
			   FieldType data) = 0;

    //! Perfom a global integration on a field for rebalance.
    virtual void integrate(std::string field_name,
			   FieldType &field_norm) = 0;

    //! Perform a rebalance on a field for global conservation.
    virtual void rebalance(std::string field_name,
			   FieldType field_norm) = 0;
};

} // end namespace coupler

#endif // coupler_Transfer_Evaluator_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Evaluator.hh
//---------------------------------------------------------------------------//
