//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   src/Transfer_Map.hh
 * \author Stuart Slattery
 * \date   Thu Oct 06 15:53:09 2011
 * \brief  Transfer_Map class definition.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef dtransfer_Transfer_Map_hh
#define dtransfer_Transfer_Map_hh

#include "Transfer_Evaluator.hh"
#include <vector>

namespace dtransfer
{

//===========================================================================//
/*!
 * \class Transfer_Map
 * \brief A basic map data structure to hold topological relationships between
 * meshes.
 */
//===========================================================================//

class Transfer_Map 
{
  public:
    //@{
    // Useful typedefs.
    typedef std::vector<int>                             Vector_Int;
    typedef const Transfer_Evaluator*                    Const_TE;
    //@}

  private:
    
    // Transfer_Evaluator pointer.
    Const_TE te;

    // Element index vector.
    Vector_Int element_index;

    // Element rank vector.
    Vector_Int element_rank;

  public:

    // Constructor.
    Transfer_Map(Const_TE te_);

    // Destructor.
    ~Transfer_Map();

    // Add an index to the map.
    void add_index(int index);

    // Add a rank to the map.
    void add_rank(int rank);

    // Return the index vector.
    const Vector_Int get_index() { return element_index; }

    // Return the rank vector.
    const Vector_Int get_rank() { return element_rank; }
};

} // end namespace dtransfer

#endif // dtransfer_Transfer_Map_hh

//---------------------------------------------------------------------------//
//              end of src/Transfer_Map.hh
//---------------------------------------------------------------------------//
