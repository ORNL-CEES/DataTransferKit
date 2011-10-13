//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   src/Data_Transfer_Manager.hh
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:44 2011
 * \brief  Data_Transfer_Manager class definiton.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef dtransfer_Data_Transfer_Manager_hh
#define dtransfer_Data_Transfer_Manager_hh

#include "Transfer_Evaluator.hh"
#include "Transfer_Map.hh"
#include "Field_DB.hh"

namespace dtransfer
{

//===========================================================================//
/*!
 * \class Data_Transfer_Manager
 * \brief The Data_Transfer_Manager manages the data transfer problem between
 * codes. 
 *
 * This initial implementation is limited to 2 codes for now but could
 * easily be expanded to many coupled codes operated by the same manager by
 * taking advantage of a common parallel topology map data
 * structure.
 */
//===========================================================================//

class Data_Transfer_Manager 
{
  public:
    
    //@{
    //! Useful typedefs.
    typedef Transfer_Evaluator::Data_Vector              Data_Vector;
    typedef Data_Vector::iterator                        Data_Iterator;
    typedef Data_Vector::const_iterator                  Const_Data_Iterator;
    //@}


  private:

    // Reference to Physics A transfer evaluator.
    Transfer_Evaluator* d_te_a;

    // Reference to Physics B transfer evaluator.
    Transfer_Evaluator* d_te_b;

    // Topology map for transfer from A to B.
    Transfer_Map* d_map_A2B;

    // Topology map for transfer from B to A.
    Transfer_Map* d_map_B2A;

    // Field database.
    Field_DB d_f_db;

  public:

    // Constructor.
    Data_Transfer_Manager(Transfer_Evaluator* TE_A_,
			  Transfer_Evaluator* TE_B_);

    // Destructor.
    ~Data_Transfer_Manager();

    // Register a field to be controlled by the manager.
    void add_field(std::string field_name);

    // Build the topology map for transfer from A to B.
    void map_A2B();

    // Build the topology map for transfer from B to A.
    void map_B2A();

    // Transfer data from A to B.
    void transfer_A2B(std::string field_name);

    // Transfer data from B to A.
    void transfer_B2A(std::string field_name);
};

} // end namespace dtransfer

#endif // dtransfer_Data_Transfer_Manager_hh

//---------------------------------------------------------------------------//
//              end of src/Data_Transfer_Manager.hh
//---------------------------------------------------------------------------//
