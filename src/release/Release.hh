//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   release/Release.hh
 * \author Thomas M. Evans
 * \date   Tue Jul 10 13:16:18 2007
 * \brief  Release function for coupler
 * \note   Copyright (C) 2007 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: Release.hh,v 1.4 2009/02/09 16:27:03 9te Exp $
//---------------------------------------------------------------------------//

#ifndef release_Release_hh
#define release_Release_hh

#include <string>

//===========================================================================//
/*!
 * \namespace coupler
 *
 * \brief Namespace that contains components used to create codes in the
 * denovo package.
 *
 */
//===========================================================================//

namespace coupler
{

const std::string release();
const std::string release_date();
const std::string build_date();
int               major_number();
int               minor_number();
int               branch_number();

}

#endif // release_Release_hh

//---------------------------------------------------------------------------//
//                        end of release/Release.hh
//---------------------------------------------------------------------------//
