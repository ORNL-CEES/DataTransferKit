//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Data_Transfer_Manager.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:44 2011
 * \brief  Data_Transfer_Manager member definitons.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Data_Transfer_Manager.hh"
#include "harness/DBC.hh"
#include <vector>
#include <algorithm>

namespace coupler
{

//---------------------------------------------------------------------------//
// Constructor.
Data_Transfer_Manager::Data_Transfer_Manager(Communicator_t comm_global,
					     Transfer_Evaluator* TE_A,
					     Transfer_Evaluator* TE_B)
    : d_comm_global(comm_global)
{
    // operate on the global communicator
    nemesis::set_internal_comm(d_comm_global);

    // Wrap the raw pointers.
    d_te_a = TE_A;
    d_te_b = TE_B;

    // Get the physics' subcommunicators.
    d_te_a->register_comm(d_comm_a);
    d_te_b->register_comm(d_comm_b);

    // Generate local to global indexers.
    if(d_te_a || d_te_b)
    {
	d_indexer_A = new LG_Indexer(d_comm_global, d_comm_a, d_te_a);
	Ensure( d_indexer_a );

	d_indexer_B = new LG_Indexer(d_comm_global, d_comm_b, d_te_b);
	Ensure( d_indexer_b );
    }

    // reset the internal communicator
    nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
// Destructor.
Data_Transfer_Manager::~Data_Transfer_Manager()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Register a field to be controlled by the manager.
void Data_Transfer_Manager::add_field(std::string field_name)
{
    // Check that the field is supported by the each code.
    Require( d_te_a->register_field(field_name) );
    Require( d_te_b->register_field(field_name) );

    // Add the new field to the database.
    d_f_db.add_field(field_name);
}

//---------------------------------------------------------------------------//
// Build the topology map for transfer from A to B for a particular field. B
// provides the target points and A is the source.
void Data_Transfer_Manager::map_A2B(std::string field_name)
{
    // Initialize map.
    d_map_A2B = new Transfer_Map();    

    // Target point coordinate vector.
    std::vector<double> points;

    // Target point handle vector.
    std::vector<Handle> handles;

    // Physics B registers its target points.
    d_te_b->register_xyz(field_name, points, handles);
    Check( points.size() % 3 == 0 );
    Check( points.size() / handles.size() == 3 );


}

//---------------------------------------------------------------------------//
// Transfer data from A to B.
void Data_Transfer_Manager::transfer_A2B(std::string field_name)
{
    // Get the field we are operating on.
    d_f_db;

    // Get the map contents.
    Transfer_Map::Vector_Int index = d_map_A2B->get_index();
    Transfer_Map::Vector_Int rank = d_map_A2B->get_rank();

    // Register the domain of A. We want to do this every time because it may
    // have changed.
    Data_Iterator domain_begin;
    Data_Iterator domain_end;
    d_te_a->register_domain(field_name, domain_begin, domain_end);
    int domain_size = domain_end - domain_begin;
    
    // Register the range of B. We want to do this every time because it may
    // have changed.
    Iterator range_begin;
    Iterator range_end;
    if (d_te_a->register_range(field_name, range_begin, range_end) )
    int range_size = range_end - range_begin;

    // Transfer A to B.
    std::vector<ValueType> data(domain_size);

    // Pull the data from A.
    std::copy(domain_begin, domain_end, data.begin());
    
    // Modify it with the map before pushing it to B.
    

    // Check the data vector size before copy.
    assert( data.size() == range_end - range_begin );
    
    // Push the data to B.
    std::copy(data.begin(), data.end(), range_begin);
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
