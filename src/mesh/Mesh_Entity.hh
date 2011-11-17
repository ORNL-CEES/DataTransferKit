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

template<class HandleType, class OrdinateType>
class Mesh_Entity 
{

};

} // end namespace mesh

#endif // mesh_Mesh_Entity_hh

//---------------------------------------------------------------------------//
//              end of mesh/Mesh_Entity.hh
//---------------------------------------------------------------------------//
