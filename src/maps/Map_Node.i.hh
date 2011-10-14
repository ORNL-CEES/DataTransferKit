//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Map_Node.i.hh
 * \author Stuart R. Slattery
 * \date   Thu May 26 15:20:22 2011
 * \brief  Map node template implementation.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Map_Node_i_hh
#define coupler_Map_Node_i_hh

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::Map_Node(
    OrdinateType        handle, 
    OrdinateType        partition, 
    const CouplingType& couple_obj)
    : d_handle(handle)
    , d_partition(partition)
    , d_couple_obj(couple_obj)
{   
    Require (d_handle >= 0);
    Require (d_partition >= 0);
}


//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Register \a data with name \a key in the node.
 */
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
void Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::add(
    const KeyType& key, 
    DataType       data)
{
    Require ( !exists(key) );

    InsertReturnType rt = d_map.insert( std::make_pair(key, data) );
    Ensure (rt.second == true);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns whether data associated with \a key exists.
 */
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
bool Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::exists(
    const KeyType& key) const
{
    return !( d_map.find(key) == d_map.end() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Retrieves data associated with \a key from the node.
 */
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
typename Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::DataType
Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::get(
    const KeyType& key) const
{
    Require ( exists(key) );

    const_iterator iter = d_map.find(key);
    Ensure ( iter != d_map.end() );

    return iter->second;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Retrieves data associated with \a key from the node.
 */
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
typename Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::DataType&
Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::get(const KeyType& key)
{
    Require ( exists(key) );

    iterator iter = d_map.find(key);
    Ensure ( iter != d_map.end() );

    return iter->second;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Removes data associated with \a key from the node.
 */
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
void Map_Node<OrdinateType_T, DataType_T, CouplingType_T>::remove(
    const KeyType& key)
{
    Require ( exists(key) );

    iterator iter = d_map.find(key);
    Ensure ( iter != d_map.end() );

    d_map.erase(iter);
}


//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Prints the data for the node with general coupling.
 */
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
void print_node(
    const Map_Node<OrdinateType_T, DataType_T, CouplingType_T> &map_node,
    std::ostream &out)
{
    out << "Not implemented for general coupling." << std::endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Prints the data for the node with Point coupling.
 */
template<class OrdinateType_T, class DataType_T>
void print_node(
    const Map_Node< OrdinateType_T, DataType_T, 
                    denovo::Point<DataType_T> > &map_node,
    std::ostream &out)
{
    typedef denovo::Point<DataType_T>                       Point_t;
    typedef Map_Node<OrdinateType_T, DataType_T, Point_t>   Map_Node_t;

    // Print the index and partition
    out << "Handle " << map_node.handle() << " on Partition " 
        << map_node.partition() << std::endl;

    // Get the point
    const Point_t &point = map_node.coupling_object();

    // Print the point
    out << std::showpoint;
    out << "   Point: (" << point.x() << ", " << point.y() << ", " << point.z() 
        << ")" << std::endl;
    out << std::noshowpoint;

    // Loop over the map and print the data
    for(typename Map_Node_t::const_iterator iter = map_node.begin(), 
                                        iter_end = map_node.end();
        iter != iter_end; ++iter)
    {
        out << "   " << iter->first << ": " << iter->second << std::endl;
    }
}


} // end namespace coupler

#endif // coupler_Map_Node_i_hh

//---------------------------------------------------------------------------//
//                 end of Map_Node.i.hh
//---------------------------------------------------------------------------//

    
