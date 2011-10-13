//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/CouplingType.hh
 * \author Gregory Davidson
 * \date   Tue Sep 20 17:20:15 2011
 * \brief  Provides types for coupling.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_CouplingType_hh
#define coupler_CouplingType_hh

namespace coupler
{

//===========================================================================//
/*!
 * \class CouplingType
 * \brief Provides the required interface for all coupling types.
 */
//===========================================================================//

template<class OrdinateType_T, class DataType_T>
class CouplingType 
{
  public:
    //@{
    //! Typedefs.
    typedef OrdinateType_T          OrdinateType;
    typedef DataType_T              DataType;
};

//===========================================================================//
/*!
 * \class PointCoupling
 * \brief Provides the required interface for point coupling.
 */
//===========================================================================//

template<class OrdinateType_T, class DataType_T, class PointType_T>
class PointCoupling : public PointType_T, CouplingType
{   };


} // end namespace coupler

#endif // coupler_CouplingType_hh

//---------------------------------------------------------------------------//
//              end of coupler/CouplingType.hh
//---------------------------------------------------------------------------//
