//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/EntityType.hh
 * \author Stuart Slattery
 * \date   Tue Oct 11 17:00:43 2011
 * \brief  Provides types for mesh entities.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef mesh_EntityType_hh
#define mesh_EntityType_hh

namespace mesh
{

//===========================================================================//
/*!
 * \class EntityType
 * \brief Provides the required interface for all entity types.
 */
//===========================================================================//

template<class OrdinateType_T, class DataType_T>
class EntityType
{
  public:
    //@{
    //! Useful typedefs.
    typedef OrdinateType_T            OrdinateType;
    typedef DataType_T                DataType;
    //@}

  public:

    //! Constructor.
    EntityType()
    { /* ... */ }

    //! Destructor.
    virtual ~EntityType()
    { /* ... */ }

    //! Return this entity's type.
    virtual void my_type(int &type) = 0;
};

} // end namespace mesh

#endif // mesh_EntityType_hh

//---------------------------------------------------------------------------//
//              end of mesh/EntityType.hh
//---------------------------------------------------------------------------//
