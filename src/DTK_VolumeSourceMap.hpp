//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file DTK_VolumeSourceMap.hpp
 * \author Stuart R. Slattery
 * \brief Volume source map declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_VOLUMESOURCEMAP_HPP
#define DTK_VOLUMESOURCEMAP_HPP

#include "DTK_GeometryTraits.hpp"
#include "DTK_GeometryManager.hpp"
#include "DTK_FieldTraits.hpp"
#include "DTK_FieldEvaluator.hpp"
#include "DTK_FieldManager.hpp"
#include "DTK_CommIndexer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Map_def.hpp>
#include <Tpetra_Directory_decl.hpp>
#include <Tpetra_Directory_def.hpp>
#include <Tpetra_Export.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class VolumeSourceMap
 * \brief A map for transfers where the volume is a source.
 *
 * In all cases supported by this map, the source consists of a group of
 * volumes distributed in parallel. That target will provide a list of
 * coordinates that data is desired for. These may be the centroids of target
 * volumes, quadrature points in a target mesh, or any other representation
 * that may be resolved by points. The apply stage will provide the source
 * with a list of local source geometries and a set of target points at which
 * to evaluate the source function. This evaluation can be interpreted in any
 * way, simply moving the volume data or doing a zero order evaluation for the
 * volume to volume case are a higher order functional evaluation for the
 * quadrature point case.
 */
//---------------------------------------------------------------------------//
template<class Geometry, class CoordinateField>
class VolumSourceMap
{
  public:

    //@{
    //! Typedefs.
    typedef Geometry                                  geometry_type;
    typedef GeometryTraits<Geometry>                  GT;
    typedef GeometryManager<Geometry>                 GeometryManagerType;
    typedef Teuchos::RCP<GeometryManagerType>         RCP_GeometryManager;
    typedef CoordinateField                           coord_field_type;
    typedef FieldTraits<CoordinateField>              CFT;
    typedef typename CFT::size_type                   CoordOrdinal;
    typedef FieldManager<CoordinateField>             CoordFieldManagerType;
    typedef Teuchos::RCP<CoordFieldManagerType>       RCP_CoordFieldManager;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;
    typedef Tpetra::Map<GlobalOrdinal>                TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    typedef Tpetra::Export<GlobalOrdinal>             ExportType;
    typedef Teuchos::RCP<ExportType>                  RCP_TpetraExport;
    //@}
};

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_VolumeSourceMap_def.hpp"

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_VOLUMESOURCEMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_VolumeSourceMap.hpp
//---------------------------------------------------------------------------//

