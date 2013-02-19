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
#include <DTK_MeshTools.hpp>
#include <DTK_SharedDomainMap.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshTraitsFieldAdapter.hpp>

#include "MoabMesh.hpp"
#include "ArrayField.hpp"
#include "PeaksEvaluator.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>


//---------------------------------------------------------------------------//
// Example driver.
//---------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
  // Typedefs.
  typedef MoabMesh::Container MeshType;
  using Teuchos::CommandLineProcessor;

  // Setup communication.
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  bool success = true;
  bool verbose = true;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  try {

    // Read command-line options
    CommandLineProcessor  clp(false); // Don't throw exceptions
    CommandLineProcessor::EParseCommandLineReturn parse_return
      = clp.parse(argc, argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
    { return parse_return; }

    // Setup source mesh. Partitioned in x.
    int mesh_dim = 2;
    MoabMesh source_mesh( comm, "tri_peaks.vtk", moab::MBTRI, 0 );
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > 
      source_blocks( 1, source_mesh.meshContainer() );
    Teuchos::RCP<DataTransferKit::MeshManager<MeshType> > 
      source_mesh_manager =
      Teuchos::rcp( new DataTransferKit::MeshManager<MeshType>(
          source_blocks, comm, mesh_dim ) );

    // Setup target coordinate field. Partitioned in y.
    MoabMesh target_mesh( comm, "quad_mesh.vtk", moab::MBQUAD, 1 );
    Teuchos::RCP<DataTransferKit::FieldManager<MeshType> > 
      target_coords_manager =
      Teuchos::rcp( new DataTransferKit::FieldManager<MeshType>(
          target_mesh.meshContainer(), comm ) );

    // Create a peaks function evaluator.
    Teuchos::RCP<DataTransferKit::FieldEvaluator<MeshType::global_ordinal_type,
      ArrayField> > 
      peaks_evaluator = Teuchos::rcp( new PeaksEvaluator( source_mesh ) );

    // Create data target.
    int num_targets = DataTransferKit::MeshTools<MeshType>::numVertices( 
      *(target_mesh.meshContainer()) );
    Teuchos::RCP<ArrayField> data_target = Teuchos::rcp(
      new ArrayField( num_targets, 1 ) );
    Teuchos::RCP<DataTransferKit::FieldManager<ArrayField> > 
      target_space_manager = 
      Teuchos::rcp( new DataTransferKit::FieldManager<ArrayField>( 
          data_target, comm ) );

    // Construct shared domain map.
    DataTransferKit::SharedDomainMap<MeshType,MeshType> 
      shared_domain_map( comm, mesh_dim );

    // Call setup. This creates the mapping.
    shared_domain_map.setup( source_mesh_manager, target_coords_manager );

    // Apply the map. This does the field evaluation and moves the data.
    shared_domain_map.apply( peaks_evaluator, target_space_manager );

    // Set the data on the target mesh and write.
    target_mesh.tag( *data_target );
    target_mesh.write( "quad_result.vtk" );
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)
  
    if (verbose) {
      if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
      else         *out << "\nOh no! At least one of the tests failed!\n";
    }

  return ( success ? 0 : 1 );

}

//---------------------------------------------------------------------------//
// end scaling_study_1.cpp
//---------------------------------------------------------------------------//

