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
    //@}

    //! Constructor.
    Transfer_Evaluator()
    { /* ... */ }

    //! Destructor.
    virtual ~Transfer_Evaluator()
    { /* ... */ }

    //! Register entities.
    virtual void register_xyz(std::vector<double> &points) = 0;

    //! Register a field associated with the entities. Return false if this
    //! field is not supported.
    virtual bool register_field(std::string field_name) = 0;

    //! Register the domain of a field.
    virtual void register_domain(std::string field_name,
	                         Const_Iterator &begin, 
				 Const_Iterator &end) = 0;

    //! Register the range of a field.
    virtual void register_range(std::string field_name,
				Iterator &begin, 
				Iterator &end) = 0;

    //! Given (x,y,z) coordinates, return the local process rank in which that
    //! point exists and the index into the local data vector that will be
    //! applied at that point. Return true if point is in the local domain,
    //! false if not. 
    virtual void find_xyz(double x, 
			  double y, 
			  double z,
			  int &rank,
			  int index,
			  bool domain) = 0;

    //! Pull data from a field.
    virtual void pull_data(std::string field_name,
			   Data_Vector &data) = 0;

    //! Push data onto a field.
    virtual void push_data(std::string field_name,
			   Data_Vector data) = 0;
};

} // end namespace dtransfer

#endif // dtransfer_Transfer_Evaluator_hh

//---------------------------------------------------------------------------//
//              end of src/Transfer_Evaluator.hh
//---------------------------------------------------------------------------//
