//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   src/Data_Transfer_Manager.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:44 2011
 * \brief  Data_Transfer_Manager member definitons.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Data_Transfer_Manager.hh"
#include <vector>
#include <algorithm>
#include <iostream>

namespace dtransfer
{

//---------------------------------------------------------------------------//
// Constructor.
Data_Transfer_Manager::Data_Transfer_Manager(Transfer_Evaluator* TE_A_,
					     Transfer_Evaluator* TE_B_)
    : d_te_a(TE_A_)
    , d_te_b(TE_B_)
{
    // Initialize maps.
    d_map_A2B = new Transfer_Map(d_te_a);
    d_map_B2A = new Transfer_Map(d_te_b);
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
    assert( d_te_a->register_field(field_name) );
    assert( d_te_a->register_field(field_name) );

    // Add the new field to the database.
    d_f_db.add_field(field_name);
}

//---------------------------------------------------------------------------//
// Build the topology map for transfer from A to B.
void Data_Transfer_Manager::map_A2B()
{
    // Point coordinate vector.
    std::vector<double> points;

    // Physics B registers its points to be mapped.
    d_te_b->register_xyz(points);

    // For every point in B, determine the topological relationship to A.
    assert( points.size() % 3 == 0 );
    std::vector<double>::const_iterator pt_it;
    int num_points = points.size() / 3;
    int rank = 0;
    int index = 0;    
    bool domain = false;
    for (pt_it = points.begin(); pt_it != points.end(); pt_it += 3)
    {
	d_te_a->find_xyz( *(pt_it),
			  *(pt_it + 1),
			  *(pt_it + 2),
			  rank,
			  index,
			  domain);

	d_map_A2B->add_rank(rank);
	d_map_A2B->add_index(index);
    }
}

//---------------------------------------------------------------------------//
// Transfer data from A to B.
void Data_Transfer_Manager::transfer_A2B(std::string field_name)
{
    // Get the map contents.
    Transfer_Map::Vector_Int index = d_map_A2B->get_index();
    Transfer_Map::Vector_Int rank = d_map_A2B->get_rank();

    // Register the domain of A. We want to do this every time because it may
    // have changed.
    Data_Iterator domain_begin;
    Data_Iterator domain_end;
    d_te_a->register_domain(field_name, domain_begin, domain_end);    
    
    // Register the range of B. We want to do this every time because it may
    // have changed.
    Iterator range_begin;
    Iterator range_end;
    d_te_a->register_range(field_name, range_begin, range_end);

    // Transfer A to B.
    ValueType data;

    // Pull the data from A.
    d_te_a->pull_data(field_name, data);
    
    // Modify it with the map before pushing it to B.
    
    // Push the data to B.
    d_te_b->push_data(field_name, data);

//    d_te_a->pull_data(index[i], data);
//    d_te_b->push_data(i, data);
}

//---------------------------------------------------------------------------//

} // end namespace dtransfer

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
