//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Base_Map.i.hh
 * \author Stuart R. Slattery
 * \date   Wed May 25 16:08:10 2011
 * \brief  Definition of class Base_Map.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Base_Map_i_hh
#define coupler_Base_Map_i_hh


namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param dest_physics  A pointer to the physics package being mapped \em to.
 */
template<class OrdinateType_T, class DataType_T>
Base_Map<OrdinateType_T, DataType_T>::Base_Map(SP_PhysicsType dest_phys)
    : b_dest_physics(dest_phys)
    , b_mapped(false)
{
    Require (dest_phys);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Print the map.
 */
template<class OrdinateType_T, class DataType_T>
void Base_Map<OrdinateType_T, DataType_T>::print(std::ostream& out) const
{
    for(typename Vec_Node::const_iterator iter = b_nodes.begin(),
                                      iter_end = b_nodes.end();
        iter != iter_end; ++iter)
    {
        iter->print(out);
    }
}

}   // end namespace coupler

#endif // coupler_Base_Map_i_hh

//---------------------------------------------------------------------------//
//                 end of Base_Map.i.hh
//---------------------------------------------------------------------------//
