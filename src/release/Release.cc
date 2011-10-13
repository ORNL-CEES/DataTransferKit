//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   release/Release.cc
 * \author Thomas M. Evans
 * \date   Tue Jul 10 13:16:18 2007
 * \brief  Release function implementation for denovo
 * \note   Copyright (C) 2007 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: Release.cc,v 1.27 2009/08/12 21:11:55 9te Exp $
//---------------------------------------------------------------------------//

#include <release/config.h>
#include "Release.hh"

namespace coupler
{

using std::string;

//---------------------------------------------------------------------------//
/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form release-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "coupler-0_0_0";
    return pkg_release;
}

//---------------------------------------------------------------------------//
/*!  
 * \return string of the release date
 */
const string release_date()
{
    string pkg_release = "3-MAR-2011";
    return pkg_release;
}

//---------------------------------------------------------------------------//
/*!
 * \return return string of build date
 */
const string build_date()
{
#ifdef BUILD_DATE
    return BUILD_DATE;
#else
    return "unknown";
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \return integer of the major release number
 */
int major_number()
{
    return 3;
}

//---------------------------------------------------------------------------//
/*!
 * \return integer of the minor release number
 */
int minor_number()
{
    return 0;
}

//---------------------------------------------------------------------------//
/*!
 * \return integer of the tag release number
 */
int branch_number()
{
    return 0;
}

}  // end of coupler

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
