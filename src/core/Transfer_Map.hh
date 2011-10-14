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

#include <vector>

namespace dtransfer
{

//===========================================================================//
/*!
 * \class Transfer_Map
 * \brief A basic map data structure to hold topological relationships between
 * meshes. The domain_index and domain_rank vector elements correspond to each
 * element in the associated range.
 */
//===========================================================================//

class Transfer_Map 
{
  public:
    //@{
    // Useful typedefs.
    typedef std::vector<int>                          Vector_Int;
    //@}

  private:

    // Domain index vector.
    Vector_Int domain_index;

    // Domain rank vector.
    Vector_Int domain_rank;

  public:

    // Constructor.
    Transfer_Map();

    // Destructor.
    ~Transfer_Map();

    // Add a rank/index pair to the map.
    void add_pair(int rank, int index);

    // Return the index vector.
    const Vector_Int get_index() { return domain_index; }

    // Return the rank vector.
    const Vector_Int get_rank() { return domain_rank; }
};

} // end namespace dtransfer

#endif // dtransfer_Transfer_Map_hh

//---------------------------------------------------------------------------//
//              end of src/Transfer_Map.hh
//---------------------------------------------------------------------------//
