//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Map_Node.hh
 * \author Stuart R. Slattery
 * \date   Thu May 26 15:20:22 2011
 * \brief  Map node template definition.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Map_Node_hh
#define coupler_Map_Node_hh

#include <string>
#include <map>
#include <ostream>

#include "mesh_type/Point.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Map_Node
 * \brief A container that holds a mapping object and all of the information
 *        associated with that object.
 *
 * \param OrdinateType_T  An ordinal type for handles and partitions.
 * \param DataType_T      The type of data associated with the mapping object.
 * \param CouplingType_T  A type indicating what coupling to use.  For instance,
 *                        for point-wise coupling, the \a CouplingType_T is a
 *                        \a Point<double>. 
 */
/*! 
 * \example coupler/test/tstMap_Node.cc
 *
 * Test of Map_Node.
 */
//===========================================================================//

template<class OrdinateType_T, class DataType_T, class CouplingType_T>
class Map_Node 
{
  public:
    //@{
    //! Useful typedefs.
    typedef OrdinateType_T                          OrdinateType;
    typedef DataType_T                              DataType;
    typedef CouplingType_T                          CouplingType;
    typedef std::string                             KeyType;
    typedef std::map<KeyType, DataType>             NodeMap;
    typedef typename NodeMap::iterator              iterator;
    typedef typename NodeMap::const_iterator        const_iterator;
    //@}

  private:
    // >>> PRIVATE TYPEDEFS
    typedef std::pair<iterator, bool>               InsertReturnType;

    // >>> PRIVATE DATA 
    // Node index.
    OrdinateType d_handle;

    // Remote application global PID that owns this node.
    OrdinateType d_partition;

    // Mapping object.
    CouplingType d_couple_obj;

    // Node map.
    NodeMap d_map;

  public:
    // Constructor
    Map_Node(OrdinateType handle, OrdinateType partition, 
             const CouplingType& couple_obj);

    // >>> ACCESSOR FUNCTIONS
    //! Get the node index.
    OrdinateType handle() const { return d_handle; }

    //! Get the parallel partition id.
    OrdinateType partition() const { return d_partition; }

    //! Get the mapping object for this node.
    const CouplingType& coupling_object() const { return d_couple_obj; }

    // Register data.
    void add(const KeyType& key, DataType data);

    // Data exists.
    bool exists(const KeyType& key) const;
     
    // Retrieve data.
    DataType get(const KeyType& key) const;

    // Retrieve data.
    DataType& get(const KeyType& key);

    // Delete data.
    void remove(const KeyType& key);

    //@{
    //! Return an iterator into the beginning of the map.
    iterator       begin()       { return d_map.begin(); }
    const_iterator begin() const { return d_map.begin(); }
    //@}

    //@{
    //! Return an iterator into the end of the map.
    iterator       end()       { return d_map.end(); }
    const_iterator end() const { return d_map.end(); }
    //@}

    // Print the data.
    void print(std::ostream& out) const { print_node(*this, out); }
};


//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
template<class OrdinateType_T, class DataType_T, class CouplingType_T>
void print_node(
    const Map_Node<OrdinateType_T, DataType_T, CouplingType_T> &map_node,
    std::ostream& out);



} // end namespace coupler

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "Map_Node.i.hh"

#endif // coupler_Map_Node_hh

//---------------------------------------------------------------------------//
//              end of coupler/Map_Node.hh
//---------------------------------------------------------------------------//
