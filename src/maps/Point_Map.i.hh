//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Point_Map.i.hh
 * \author Stuart R. Slattery
 * \date   Wed May 25 16:51:10 2011
 * \brief  Point_Map member definitions
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include <set>
#include <functional>

#include "harness/DBC.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor taking a partitioner.
 */
template<class OrdinateType_T, class DataType_T>
Point_Map<OrdinateType_T, DataType_T>::Point_Map(SP_Partitioner part)
    : Base(part)
{ 
    Require (b_partitioner);
}


//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*! 
 * \brief Add a point to the map.
 */
template<class OrdinateType_T, class DataType_T>
void Point_Map<OrdinateType_T, DataType_T>::add_point(const Point_t &point)
{
    Require (b_partitioner);

    // Get the next index number
    OrdinateType index = b_nodes.size();

    // ask neutronics partitioner for the point's partition
    OrdinateType partition = b_partitioner->point_query( point.coords() );

    // add the new map node to the node vector
    b_nodes.push_back( Map_Node_t(index, partition, point) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Insert a power value into the data field.
 */
template<class OrdinateType_T, class DataType_T>
void Point_Map<OrdinateType_T, DataType_T>::insert_data(
    OrdinateType index, const KeyType &key, DataType data)
{   
    Require (b_mapped);

    iterator iter = find_map_node(index);
    Check (iter != b_nodes.end());

    // Insert the power
    iter->add(key, data);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Retrieve data from a node.
 */
template<class OrdinateType_T, class DataType_T>
typename Point_Map<OrdinateType_T, DataType_T>::DataType
Point_Map<OrdinateType_T, DataType_T>::get_data(OrdinateType index,
                                                const KeyType &key) const
{
    Require (b_mapped);

    // Find the node
    const_iterator iter = find_map_node(index);
    Check ( iter != b_nodes.end() );
    
    // Return the data
    return iter->get(key);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Retrieve data from a node.
 */
template<class OrdinateType_T, class DataType_T>
typename Point_Map<OrdinateType_T, DataType_T>::DataType&
Point_Map<OrdinateType_T, DataType_T>::get_data(OrdinateType index,
                                                const KeyType &key)
{
    Require (b_mapped);

    // Find the node
    iterator iter = find_map_node(index);
    Check ( iter != b_nodes.end() );
    
    // Return the data
    return iter->get(key);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns whether a given \a key exists for node \a index.
 */
template<class OrdinateType_T, class DataType_T>
bool Point_Map<OrdinateType_T, DataType_T>::data_exists(
    OrdinateType index, const KeyType &key) const
{
    Require (b_mapped);

    // Find the node
    const_iterator iter = find_map_node(index);

    if( iter == b_nodes.end() )
        return false;
    else
        return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Register the mapping as complete.
 */
template<class OrdinateType_T, class DataType_T>
void Point_Map<OrdinateType_T, DataType_T>::complete()
{
    // do a binary sort on the nodes by partition id
    sort_nodes_by_partition();

    // Make a vector of which partitions are in the nodes
    b_partitions.push_back( b_nodes.front().partition() );
    for(const_iterator iter = b_nodes.begin()+1, iter_end = b_nodes.end();
        iter != iter_end; ++iter)
    {
        if( iter->partition() != b_partitions.back() )
            b_partitions.push_back( iter->partition() );
    }

    // set up the partition_pts member vector
    b_partition_pts.resize(b_partitions.size(), 0);

    // calculate the number of local points on each remote application
    // partition
    int n = 0;
    for (iterator iter = b_nodes.begin(), iter_end = b_nodes.end(); 
         iter != iter_end; ++iter)
    {
        if ( iter->partition() == b_partitions[n] )
        {
            ++b_partition_pts[n];
        }
        else
        {
            ++n;
            ++b_partition_pts[n];
        }
    }

    // signal that mapping is complete
    b_mapped = true;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get a map node with a specific index.
 */
template<class OrdinateType_T, class DataType_T>
typename Point_Map<OrdinateType_T, DataType_T>::Map_Node_t& 
Point_Map<OrdinateType_T, DataType_T>::get_node(OrdinateType index)
{
    // Get the map node with the given index
    iterator iter = find_map_node(index);
    Check ( iter != b_nodes.end() );

    return *iter;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a map node with a specific index.
 */
template<class OrdinateType_T, class DataType_T>
const typename Point_Map<OrdinateType_T, DataType_T>::Map_Node_t&
Point_Map<OrdinateType_T, DataType_T>::get_node(OrdinateType index) const
{
    // Get the map node with the given index
    const_iterator iter = find_map_node(index);
    Check ( iter != b_nodes.end() );

    return *iter;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sort the node vector by partition.
 */
template<class OrdinateType_T, class DataType_T>
void Point_Map<OrdinateType_T, DataType_T>::sort_nodes_by_partition()
{
    std::sort(b_nodes.begin(), b_nodes.end(), 
              &Point_Map<OrdinateType, DataType>::partition_less_than);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sort the node vector by handle.
 */
template<class OrdinateType_T, class DataType_T>
void Point_Map<OrdinateType_T, DataType_T>::sort_nodes_by_handle()
{
    std::sort(b_nodes.begin(), b_nodes.end(), 
              &Point_Map<OrdinateType, DataType>::handle_less_than);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the range of nodes with the given partition.
 */
template<class OrdinateType_T, class DataType_T>
typename Point_Map<OrdinateType_T, DataType_T>::iterator_pair
Point_Map<OrdinateType_T, DataType_T>::find_node_range(OrdinateType partition)
{
    typedef NodeComparator<OrdinateType_T, DataType_T>  NodeComparator_t;

    // Get the beginning of the range
    iterator begin = 
        std::find_if(b_nodes.begin(), b_nodes.end(), 
                     NodeComparator_t(NodeComparator_t::PARTITION_MODE, 
                                      partition) );

    if( begin == b_nodes.end() )
    {
        return std::make_pair(begin, begin);
    }
    else
    {
        // Get the end of the range
        iterator end = 
            std::find_if(begin, b_nodes.end(), 
                         NodeComparator_t(NodeComparator_t::PARTITION_MODE, 
                                          partition) );

        return std::make_pair(begin, end);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the range of nodes with the given partition.
 */
template<class OrdinateType_T, class DataType_T>
typename Point_Map<OrdinateType_T, DataType_T>::const_iterator_pair
Point_Map<OrdinateType_T, DataType_T>::find_node_range(
    OrdinateType partition) const
{
    typedef NodeComparator<OrdinateType_T, DataType_T>  NodeComparator_t;

    // Get the beginning of the range
    const_iterator begin = 
        std::find_if(b_nodes.begin(), b_nodes.end(), 
                     NodeComparator_t(NodeComparator_t::PARTITION_MODE, 
                                      partition) );

    if( begin == b_nodes.end() )
    {
        return std::make_pair(begin, begin);
    }
    else
    {
        // Get the end of the range
        const_iterator end = 
            std::find_if(begin, b_nodes.end(), 
                         NodeComparator_t(NodeComparator_t::PARTITION_MODE, 
                                          partition) );

        return std::make_pair(begin, end);
    }
}


//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Find a node with a particular index.
 */
template<class OrdinateType_T, class DataType_T>
typename Point_Map<OrdinateType_T, DataType_T>::iterator
Point_Map<OrdinateType_T, DataType_T>::find_map_node(OrdinateType index)
{
    typedef NodeComparator<OrdinateType_T, DataType_T>      NodeComparator_t;

    return std::find_if(b_nodes.begin(), b_nodes.end(), 
                        NodeComparator_t(NodeComparator_t::INDEX_MODE, index) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a node with a particular index.
 */
template<class OrdinateType_T, class DataType_T>
typename Point_Map<OrdinateType_T, DataType_T>::const_iterator
Point_Map<OrdinateType_T, DataType_T>::find_map_node(OrdinateType index) const
{
    typedef NodeComparator<OrdinateType_T, DataType_T>      NodeComparator_t;

    return std::find_if(b_nodes.begin(), b_nodes.end(), 
                        NodeComparator_t(NodeComparator_t::INDEX_MODE, index) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the less-than operator on the partitions of two nodes.
 */
template<class OrdinateType_T, class DataType_T>
bool Point_Map<OrdinateType_T, DataType_T>::partition_less_than(
    const Map_Node_t &n1, const Map_Node_t &n2)
{
    return n1.partition() < n2.partition();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the less-than operator on the handles of two nodes.
 */
template<class OrdinateType_T, class DataType_T>
bool Point_Map<OrdinateType_T, DataType_T>::handle_less_than(
    const Map_Node_t &n1, const Map_Node_t &n2)
{
    return n1.point().handle() < n2.point().handle();
}


}   // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Point_Map.cc
//---------------------------------------------------------------------------//
