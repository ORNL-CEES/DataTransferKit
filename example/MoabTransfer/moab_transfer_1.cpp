//---------------------------------------------------------------------------//
/*!
 * \file moab_transfer_1.cpp
 * \author Stuart R. Slattery
 * \brief Moab mesh transfer example.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include <DTK_MeshManager.hpp>
#include <DTK_TransferOperator.hpp>
#include <DTK_SharedDomainMap.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTraitsFieldAdapter.hpp>

#include "MoabMesh.hpp"
#include "ArrayField.hpp"
#include "PeaksEvaluator.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>

//---------------------------------------------------------------------------//
// Example driver.
//---------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Setup source mesh.
    typedef MeshTraits<MoabMesh> MT;
    int mesh_dim = 2;
    MoabMesh source_mesh( comm, "tri_peaks.vtk", moab::MBTRI, 0 );
    Teuchos::ArrayRCP<MoabMesh> source_blocks( 1, source_mesh );
    Teuchos::RCP< DataTransferKit::MeshManager<MoabMesh> > source_mesh_manager =
	Teuchos::rcp( new DataTransferKit::MeshManager<MoabMesh>(
			  source_blocks, comm, mesh_dim ) );

    // Setup target coordinate field.
    MoabMesh target_mesh( comm, "quad_mesh.vtk", moab::MBQUAD, 1 );
    Teuchos::RCP< DataTransferKit::FieldManager<MT> > target_coords_manager =
	Teuchos::rcp( new DataTransferKit::FieldManager<MT>(
			  target_mesh, comm ) );

    // Create a peaks function evaluator.
    Teuchos::RCP< DataTransferKit::FieldEvaluator<MoabMesh,ArrayField> > 
	peaks_evaluator = Teuchos::rcp( new PeaksEvaluator( source_mesh ) );

    // Create data target.
    ArrayField data_target( target_coords.size() / mesh_dim, 1 );

    // Setup shared domain mapping.
    typedef DataTransferKit::SharedDomainMap<MoabMesh,MoabMesh> MapType;
    Teuchos::RCP<MapType> shared_domain_map = 
    	Teuchos::rcp( new MapType( comm ) );

    // Create the transfer operator.
    DataTransferKit::TransferOperator<MapType> 
	transfer_operator( shared_domain_map );

    // Setup the transfer operator ( this creates the mapping ).
    transfer_operator.setup( source_mesh_manager, target_coords_manager );

    // Apply the transfer operator ( this does the field evaluation and moves
    // the data ).
    transfer_operator.apply( peaks_evaluator, target_coords );

    // Set the data on the target mesh and write.
    target_mesh.tag( data_target );
    target_mesh.write( "quad_result.vtk" );

    return 0;
}

//---------------------------------------------------------------------------//
// end scaling_study_1.cpp
//---------------------------------------------------------------------------//

