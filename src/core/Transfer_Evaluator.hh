//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   src/Transfer_Evaluator.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:59 2011
 * \brief  Transfer_Evaluator class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef dtransfer_Transfer_Evaluator_hh
#define dtransfer_Transfer_Evaluator_hh

#include "comm/global.hh"
#include <vector>
#include <string>

namespace dtransfer
{

//===========================================================================//
/*!
 * \class Transfer_Evaluator
 * \brief Base class definition of the transfer evaluator interface.
 */
//===========================================================================//
template<class FieldType_T>
class Transfer_Evaluator 
{
  public:

    //@{
    //! Useful typedefs.
    typedef FieldType_T                                   FieldType;
    typedef typename FieldType::value_type                ValueType;
    typedef typename ValueType::iterator                  Iterator;
    typedef typename ValueType::const_iterator            Const_Iterator;
    typedef nemesis::Communicator_t                       Communicator_t;
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

    //! Register entities with a field.
    virtual void register_xyz(std::string field_name,
			      std::vector<double> &points) = 0;

    //! Register the domain of a field. Return false if domain not supported
    //! for the field.
    virtual void register_domain(std::string field_name,
	                         Const_Iterator &begin, 
				 Const_Iterator &end) = 0;

    //! Register the range of a field. Return false if range not supported for
    //! the field.
    virtual void register_range(std::string field_name,
				Iterator &begin, 
				Iterator &end) = 0;

    //! Given (x,y,z) coordinates, return the process rank on which that
    //! point exists and the iterator into the local data vector that will be
    //! applied at that point. Return true if point is in the local domain,
    //! false if not. 
    virtual bool find_xyz(double x, 
			  double y, 
			  double z,
			  int &rank,
			  Const_Iterator &value_iterator) = 0;

    //! Perfom a global integration on a field for rebalance.
    virtual void integrate_field(std::string field_name,
				 ValueType &field_norm) = 0;

    //! Perform a rebalance on a field for global/local conservation.
    virtual void rebalance(std::string field_name,
			   ValueType field_norm) = 0;
};

} // end namespace dtransfer

#endif // dtransfer_Transfer_Evaluator_hh

//---------------------------------------------------------------------------//
//              end of src/Transfer_Evaluator.hh
//---------------------------------------------------------------------------//
