//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Mesh_Entity.hh
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:59:04 2011
 * \brief  Base interface for mesh entities.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef mesh_Mesh_Entity_hh
#define mesh_Mesh_Entity_hh

namespace mesh
{

//===========================================================================//
/*!
 * \class Mesh_Entity
 * \brief Base interface for mesh entities.
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at http://www.doxygen.org.
 *
 * \sa Mesh_Entity.cc for detailed descriptions.
 *
 * \par Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*! 
 * \example mesh/test/tstMesh_Entity.cc
 *
 * Test of Mesh_Entity.
 */
//===========================================================================//

template<class HandleType_T, class OrdinateType_T>
class Mesh_Entity 
{

  public:

    //@{
    //! Useful typedefs.
    typedef HandleType_T            HandleType;
    typedef OrdinateType_T          OrdinateType;
    //@}

  protected:

    // Entity handle.
    HandleType b_handle.
    
  public:

    // Constructor.
    Mesh_Entity()
    { /* ... */ }

    // Destructor.
    virtual ~Mesh_Entity()
    { /* ... */ }

    // Return the handle of this entity.
    HandleType handle() { return b_handle; }

    // Return the dimension of this entity.
    virtual int dimension()
    { return -1; }

};

} // end namespace mesh

#endif // mesh_Mesh_Entity_hh

//---------------------------------------------------------------------------//
//              end of mesh/Mesh_Entity.hh
//---------------------------------------------------------------------------//
