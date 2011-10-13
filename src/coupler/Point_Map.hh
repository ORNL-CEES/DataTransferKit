//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Point_Map.hh
 * \author Stuart R. Slattery
 * \date   Wed May 25 16:51:10 2011
 * \brief  Point_Map class definition
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Point_Map_hh
#define coupler_Point_Map_hh

#include <vector>
#include <algorithm>
#include <functional>

#include "mesh_type/Point.hh"
#include "Map_Node.hh"
#include "Base_Map.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Point_Map
 * \brief Realization of Base_Map class for field point coupling
 *
 * \sa Point_Map.cc for detailed descriptions.
 *
 * \sa tstPoint_Map.cc for usage examples.
 */
/*! 
 * \example coupler/test/tstPoint_Map.cc
 *
 * Test of Point_Map.
 */
//===========================================================================//

template<class DataType_T, class PhysicsType_T, class PointType_T>
class Point_Map : public Base_Map<OrdinateType_T, DataType_T,
                                  PointCoupling<>
{
  public:
    //@{
    //! Useful typedefs.
    typedef Base_Map<OrdinateType_T, DataType_T>        Base;
    typedef OrdinateType_T                              OrdinateType;
    typedef DataType_T                                  DataType;
    typedef denovo::Point<OrdinateType_T, DataType_T>   Point_t;
    typedef typename Base::SP_Partitioner               SP_Partitioner;
    typedef typename Base::Map_Node_t                   Map_Node_t;
    typedef typename Base::KeyType                      KeyType;
    typedef typename Base::Vec_Node                     Vec_Node;
    typedef typename Base::Vec_Ord                      Vec_Ord;
    typedef typename Base::iterator                     iterator;
    typedef typename Base::const_iterator               const_iterator;
    typedef typename Base::iterator_pair                iterator_pair;
    typedef typename Base::const_iterator_pair          const_iterator_pair;
    //@}

    //@{
    //! Expose base class members.
    using Base::b_mapped;
    using Base::b_partitioner;
    using Base::b_nodes;
    using Base::b_partitions;
    using Base::b_partition_pts;
    //@}

    // Constructor
    explicit Point_Map(SP_Partitioner part);

    // Destructor
    ~Point_Map() {  }
    
    // Add a point to the map
    void add_point(const Point_t &point);

    // Insert data into a map node
    virtual void insert_data(OrdinateType index, const KeyType &key,
                             DataType data);
    
    // Retrieve data from a node.
    virtual DataType  get_data(OrdinateType index, const KeyType &key) const;
    virtual DataType& get_data(OrdinateType index, const KeyType &key);

    // Returns whether a given \a key exists for node \a index.
    virtual bool data_exists(OrdinateType index, const KeyType &key) const;

    // Complete the mapping
    virtual void complete();

    // Get a node with a specific point handle
    virtual Map_Node_t& get_node(OrdinateType index);
    virtual const Map_Node_t& get_node(OrdinateType index) const;

    // Sort the map node vector by node partition.
    virtual void sort_nodes_by_partition();
    // Sort the map node vector by node point handle.
    virtual void sort_nodes_by_handle();
    // Find the range of nodes with the given partition
    virtual iterator_pair find_node_range(OrdinateType partition);
    virtual const_iterator_pair find_node_range(OrdinateType partition) const;

  private:
    // >>> IMPLEMENTATION
    // Find a node with a particular index.
    iterator       find_map_node(OrdinateType index);
    const_iterator find_map_node(OrdinateType index) const;

    // Less-than comparator function for standard sort of nodes by partition. 
    static bool partition_less_than(const Map_Node_t &n1, const Map_Node_t &n2);
    // Comparator function for standard sort on nodes by node point handle.
    static bool handle_less_than(const Map_Node_t &n1, const Map_Node_t &n2);
};


//===========================================================================//
/*!
 * \class NodeComparator
 * \brief Compares a node to the provided data.
 *
 * \param INDEX_MODE  In this mode, a node's index is compared to the provided
 *                    data. 
 * \param PARTITION_MODE In this mode, a node's partition is compared to the
 *                       provided data.
 */
//===========================================================================//
 /*
template<typename OrdinateType_T, class DataType_T>
class NodeComparator
{
  public:
    //@{
    //! Typedefs.
    typedef OrdinateType_T                          OrdinateType;
    typedef Point_Map<OrdinateType_T, DataType_T>   Point_Map_t;
    typedef typename Point_Map_t::Map_Node_t        Map_Node_t;
    //@}

    //! Enumeration.
    enum Mode { INDEX_MODE, PARTITION_MODE };

  private:
    // Stores the index being compared to.
    OrdinateType d_comparison_data;
    // Stores the mode of comparison being performed.
    Mode d_mode;

  public:
    //! Constructor.
    NodeComparator(Mode mode, OrdinateType comparison_data)
        : d_comparison_data(comparison_data)
        , d_mode(mode)
    {   }

    //! Comparator operator
    bool operator()(const Map_Node_t &node) const
    {
        if(d_mode == INDEX_MODE)
            return node.index() == d_comparison_data;
        else
            return node.partition() == d_comparison_data;
    }
};
*/

} // end namespace coupler

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

//#include "Point_Map.i.hh"

#endif // coupler_Point_Map_hh

//---------------------------------------------------------------------------//
//              end of coupler/Point_Map.hh
//---------------------------------------------------------------------------//
