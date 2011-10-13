//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Base_Map.hh
 * \author Stuart R. Slattery
 * \date   Wed May 25 16:08:10 2011
 * \brief  Definition of class Base_Map.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Base_Map_hh
#define coupler_Base_Map_hh

#include <vector>

#include "Map_Node.hh"
#include "kba_mesh/Partitioner.hh"
#include "utils/SP.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Base_Map
 * \brief Abstract base class for neutronics coupling maps.
 */
//===========================================================================//
template<class DataType_T, class PhysicsType_T, class CouplingType_T>
class Base_Map 
{
  public:
    //@{
    //! Useful typedefs.
    typedef int                                             OrdinateType;
    typedef DataType_T                                      DataType;
    typedef PhysicsType_T                                   PhysicsType;
    typedef CouplingType_T                                  CouplingType;
    typedef denovo::SP<PhysicsType>                         SP_PhysicsType;
    typedef Map_Node<OrdinateType, DataType, CouplingType>  Map_Node_t;
    typedef typename Map_Node_t::KeyType                    KeyType;
    typedef std::vector<Map_Node_t>                         Vec_Node;
    typedef std::vector<OrdinateType_T>                     Vec_Ord;
    typedef typename Vec_Node::iterator                     iterator;
    typedef typename Vec_Node::const_iterator               const_iterator;
    typedef std::pair<iterator, iterator>                   iterator_pair;
    typedef std::pair<const_iterator, const_iterator>       const_iterator_pair;
    //@}

  protected:
    // Stores a pointer to the physics being mapped upon.
    SP_PhysicsType b_dest_physics;
    
    // Switch signaling the mapping has been completed and data
    // operations can be performed.
    bool b_mapped;

    // Local data field map nodes.
    Vec_Node b_nodes;

    // Unique set of partitions for the remote application.
    Vec_Ord b_partitions;

    // The number of local points on each remote application partition.
    Vec_Ord b_partition_pts;

  public:
    // Constructor.
    explicit Base_Map(SP_PhysicsType dest_physics);
    
    //! Virtual destructor.
    virtual ~Base_Map() 
    {   }

    // >>> INTERFACE
    //! Insert data into map node \a index with \a key.
    virtual void insert_data(OrdinateType handle, const KeyType &key, 
                             DataType data) = 0;

    //! Retrieve data from map node \a handle with \a key.
    virtual DataType get_data(OrdinateType handle, const KeyType &key) const = 0;

    //! Retrieve data from map node \a handle with \a key.
    virtual DataType& get_data(OrdinateType handle, const KeyType &key) = 0;

    //! Returns whether a given \a key exists for node \a handle.
    virtual bool data_exists(OrdinateType handle, const KeyType &key) const = 0;

    //! Complete the mapping.
    virtual void complete() = 0;

    //@{
    //! Return a node with a specific \a handle.
    virtual Map_Node_t& get_node(OrdinateType handle) = 0;
    virtual const Map_Node_t& get_node(OrdinateType handle) const = 0;
    //@}

    //! Return the number of nodes in the map.
    typename Vec_Node::size_type num_nodes() const { return b_nodes.size(); }

    //@{
    //! Return the vector of nodes.
    const Vec_Node& nodes() const { return b_nodes; }
    Vec_Node& nodes() { return b_nodes; }
    //@}

    //! Get the unique set of remote application partition ids.
    const Vec_Ord& partitions() const { return b_partitions; }

    //! Get the number of local points on each remote application partition.
    const Vec_Ord& partition_pts() const { return b_partition_pts; }

    //! Get the state of the map.
    bool status() const { return b_mapped; }

    // Print the map
    void print(std::ostream& out) const;

    // >>> SORTING AND FINDING
    //! Sort the map by partition.
    virtual void sort_nodes_by_partition() = 0;
    //! Sort the map by handle.
    virtual void sort_nodes_by_handle() = 0;
    //! Find the range of nodes that belong to a partition.
    virtual iterator_pair find_node_range(OrdinateType partition) = 0;
    //! Find the range of nodes that belong to a partition.
    virtual const_iterator_pair find_node_range(
        OrdinateType partition) const = 0;

    // >>> ITERATORS
    //@{
    //! Return an iterator to the beginning of the nodes.
    iterator       begin()       { return b_nodes.begin(); }
    const_iterator begin() const { return b_nodes.begin(); }
    //@}

    //@{
    //! Return an iterator to the end of the nodes.
    iterator       end()         { return b_nodes.end(); }
    const_iterator end() const   { return b_nodes.end(); }
    //@}
};


} // end namespace coupler

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "Base_Map.i.hh"

#endif // coupler_Base_Map_hh

//---------------------------------------------------------------------------//
//              end of coupler/Base_Map.hh
//---------------------------------------------------------------------------//
