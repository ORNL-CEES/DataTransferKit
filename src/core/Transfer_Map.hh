//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Map.hh
 * \author Stuart Slattery
 * \date   Thu Oct 06 15:53:09 2011
 * \brief  Transfer_Map class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Transfer_Map_hh
#define coupler_Transfer_Map_hh

#include <map>
#include <set>
#include <string>

namespace coupler
{

//===========================================================================//
/*!
 * \class Transfer_Map
 * \brief A basic map data structure to hold topological relationships between
 * parallel meshes.
 *
 * For one way transfer, the source application must know which of its local
 * mesh entities, identified here by a handle, must map onto a mesh entity in
 * the target application and on what parallel process rank of that target
 * application that entity exists. The target application must know the
 * inverse of this; for its local target entities, it needs to know which
 * source processes will be sending data to map onto those entities, and which
 * pieces of data correspond to those entities.
 *
 * In the nomenclature used here, the domain of a target rank consists of the
 * mesh entity handles in the source domain that belong to that target
 * rank. The range of a source rank consists of the mesh entity handles in the
 * target domain that belong to that source rank.
 */
//===========================================================================//

class Transfer_Map 
{
  public:

    //@{
    //! Useful typedefs.
    typedef int                                          HandleType;
    typedef int                                          OrdinateType;
    typedef std::pair<OrdinateType,HandleType>           Map_Element;
    typedef std::multimap<OrdinateType,HandleType>       Map;
    typedef typename Map::const_iterator                 Map_Iterator;
    typedef std::pair<Map_Iterator,Map_Iterator>         Map_Pair;
    typedef std::set<OrdinateType>                       Set;
    typedef typename Set::const_iterator                 Set_Iterator;
    typedef std::pair<Set_Iterator,Set_Iterator>         Set_Pair;
    //@}

  private:

    // Map for the source application. 
    // Key: target rank, Value: source handle. 
    Map source_map;

    // Map for the target application. 
    // Key: source rank, Value: target handle.
    Map target_map;

    // Set of unique source ranks.
    Set source_set;

    // Set of unique target ranks.
    Set target_set;

  public:

    //! Constructor.
    inline Transfer_Map();

    //! Destructor.
    inline ~Transfer_Map();

    //! Add a pair to the source map. The source handle corresponding to an
    //! entity in the local domain correlates to the range owned by the
    //! target rank.
    inline void add_domain_pair(OrdinateType target_rank, 
				HandleType source_handle);

    //! Add a pair to the target map. The target handle corresponding to an
    //! entity in the local range correlates to the domain owned by the
    //! source rank.
    inline void add_range_pair(OrdinateType source_rank, 
			       HandleType target_handle);

    //! Get the number of source handles with a specific target rank. This is
    //! the size of the local domain correlating to the target rank.
    inline int domain_size(OrdinateType target_rank);
    
    //! Get the number of target handles with a specific source rank. This is
    //! the size of the local range correlating to the source rank.
    inline int range_size(OrdinateType source_rank);

    //! Get the iterator pair for the source domain of a target rank. This
    //! correlates to the local source handles that exist in the range of the
    //! target rank.
    inline Map_Pair domain(OrdinateType target_rank);

    //! Get the iterator pair for the target range of a source rank. This
    //! correlates to the local target handles that exist in the domain of the
    //! source rank.
    inline Map_Pair range(OrdinateType source_rank);

    //! Return a const_iterator pair to the beginning and end of the source
    //! rank set. 
    inline Set_Pair source_set();

    //! Return a const_iterator pair to the beginning and end of the target
    //! rank set. 
    inline Set_Pair target_set();
;

} // end namespace coupler

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Transfer_Map.i.hh"

#endif // coupler_Transfer_Map_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Map.hh
//---------------------------------------------------------------------------//
