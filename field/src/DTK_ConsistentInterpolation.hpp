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
 * \file DTK_ConsistentInterpolation.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation mapping declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATION_HPP
#define DTK_CONSISTENTINTERPOLATION_HPP

#include "DTK_FieldTraits.hpp"
#include "DTK_FieldEvaluator.hpp"
#include <DTK_MeshTraits.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Export.hpp>

namespace DataTransferKit
{

template<class Mesh, class CoordinateField>
class ConsistentInterpolation
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                      mesh_type;
    typedef MeshTraits<Mesh>                          MT;
    typedef typename MT::global_ordinal_type          global_ordinal_type;
    typedef CoordinateField                           coord_field_type;
    typedef FieldTraits<CoordinateField>              CFT;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;
    typedef Tpetra::Map<global_ordinal_type>          TpetraMap;
    typedef Teuchos::RCP<const TpetraMap>             RCP_TpetraMap;
    typedef Tpetra::Export<global_ordinal_type>       ExportType;
    typedef Teuchos::RCP<ExportType>                  RCP_Export;
    //!@}

    // Constructor.
    ConsistentInterpolation( const RCP_Comm& comm );

    // Destructor.
    ~ConsistentInterpolation();

    // Setup for interpolation.
    void setup( const Mesh& mesh, const CoordinateField& coordinate_field );

    // Apply the interpolation.
    template<class SourceField, class TargetField>
    void apply( const Teuchos::RCP< 
		    FieldEvaluator<Mesh,SourceField> >& source_evaluator,
		TargetField& target_space );

  private:

    // Build the bounding box for the rendezvous decomposition.
    BoundingBox buildRendezvousBox( const Mesh& mesh, 
				    const BoundingBox& mesh_box,
				    const CoordinateField& coordinate_field, 
				    const BoundingBox& coord_box );

    // Compute globally unique ordinals for the points in the coordinate
    // field.
    Teuchos::Array<global_ordinal_type> computePointOrdinals(
	const CoordinateField& coordinate_field );

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Export field map.
    RCP_TpetraMap d_export_map;

    // Import field map.
    RCP_TpetraMap d_import_map;

    // Field data exporter.
    RCP_Export d_data_export;

    // Local source elements.
    Teuchos::Array<global_ordinal_type> d_source_elements;

    // Local target coords.
    Teuchos::Array<double> d_target_coords;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_ConsistentInterpolation_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_CONSISTENTINTERPOLATION_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolation.hpp
//---------------------------------------------------------------------------//

