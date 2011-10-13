//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Coupler.i.hh
 * \author Stuart R. Slattery
 * \date   Fri Jun 10 08:56:27 2011
 * \brief  Member definitions of class Coupler.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.i.hh,v 1.4 2008/01/04 22:50:12 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Coupler_i_hh
#define coupler_Coupler_i_hh

namespace coupler
{


//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Point handle comparator function.
 */
template<class OrdinateType_T, class DataType_T>
bool Coupler<OrdinateType_T, DataType_T>::handle_compare(Point_t p1, Point_t p2)
{
    return p1.get_handle() < p2.get_handle();
}

} // end namespace coupler

#endif // coupler_Coupler_i_hh

//---------------------------------------------------------------------------//
//              end of coupler/Coupler.i.hh
//---------------------------------------------------------------------------//
