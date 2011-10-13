//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/STAR_Data.hh
 * \author Stuart R. Slattery
 * \date   Tue Jun 07 11:14:56 2011
 * \brief  Star-CCM+ Data container class definition
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Star_Data_hh
#define coupler_Star_Data_hh

#include "mesh_type/Point.hh"
#include "harness/DBC.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class STAR_Data
 * \brief Holds STAR-CCM+ data for coupling.
 */
/*! 
 * \example coupler/test/tstStar_Data.cc
 *
 * Test of Star_Data.
 */
//===========================================================================//

template<class DataType>
class STAR_Data 
{
  public:
    //@{ 
    //! Useful typedefs.
    typedef DataType                       value_type;
    typedef denovo::Point<value_type>      Point_vt;
    //@}

    //! Constructor.
    STAR_Data(Point_vt pt, value_type vol)
        : d_point(pt)
        , d_volume(vol)
    {  }

    //@{
    //! Get the point.
    Point_vt  point() const { return d_point; }
    Point_vt& point()       { return d_point; }
    //@}

    //@{
    //! Get the volume.
    value_type  volume() const { return d_volume; }
    value_type& volume()       { return d_volume; }
    //@}

  private:
    // Star centroid point
    Point_vt d_point;

    // Star cell volume
    value_type d_volume;
};

} // end namespace coupler

#endif // coupler_Star_Data_hh

//---------------------------------------------------------------------------//
//              end of coupler/Star_Data.hh
//---------------------------------------------------------------------------//
