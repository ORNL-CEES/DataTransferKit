//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Point.hh
 * \author Stuart Slattery
 * \date   Wed Oct 12 12:00:03 2011
 * \brief  Provides the required interface for points.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef mesh_Point_hh
#define mesh_Point_hh

#include "EntityType.hh"

namespace mesh
{

//===========================================================================//
/*!
 * \class Point
 * \brief Provides the required interface for points.
 */
//===========================================================================//

template<class HandleType_T, class DataType_T>
class Point : public EntityType
{
  public:
    //@{
    //! Useful Typedefs.
    typedef HandleType_T                     HandleType;
    typedef DataType_T                       DataType;
    //@}

  private:

  public:

    //! Constructor.
    Point(HandleType handle, 
	  DataType x, 
	  DataType y, 
	  DataType z);

    //! Destructor.
    ~Point();

    //! Return this entity's type.
    void my_type(int &type);
};

} // end namespace mesh

#endif // mesh_Point_hh

//---------------------------------------------------------------------------//
//              end of mesh/Point.hh
//---------------------------------------------------------------------------//
