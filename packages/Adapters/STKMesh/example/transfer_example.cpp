//---------------------------------------------------------------------------//
/*!
 * \file   transfer_example.cpp
 * \author Stuart Slattery
 * \brief  File-based data transfer driver.
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

#include <Teuchos_GlobalMPISession.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_Basis.hpp>
#include <Intrepid_CellTools.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_io/IossBridge.hpp>

#include <Ionit_Initializer.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_Region.h>

#include <boost/tr1/unordered_map.hpp>

#include "Topaz_BasisFactory.hpp"
#include "Topaz_Tolerances.hpp"
#include "Topaz_TransferOperator.hpp"
#include "Topaz_FieldEvaluator.hpp"
#include "Topaz_IntrepidSideCell.hpp"
#include "Topaz_SurfaceMesh.hpp"
#include "Topaz_STKMeshInterface.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
    return Teuchos::DefaultComm<Ordinal>::getComm();
}

//---------------------------------------------------------------------------//
// Register the element blocks in an exodus file.
void registerElementBlocks( 
    const Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData>& meta_data,
    stk_classic::io::MeshData& mesh_data,
    stk_classic::mesh::PartVector& element_block_parts )
{
    element_block_parts.clear();
    const Ioss::ElementBlockContainer& elem_blocks = 
	mesh_data.m_input_region->get_element_blocks();

    Ioss::ElementBlockContainer::const_iterator itr;
    for ( itr = elem_blocks.begin(); itr != elem_blocks.end(); ++itr) 
    {
	Ioss::GroupingEntity* entity = *itr;
	const std::string& name = entity->name();
	element_block_parts.push_back( meta_data->get_part(name) );
   }
}

//---------------------------------------------------------------------------//
// Register the sidesets in an exodus file.
void registerSidesets(
    const Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData>& meta_data, 
    stk_classic::io::MeshData& mesh_data,
    stk_classic::mesh::PartVector& sideset_parts )
{
    sideset_parts.clear();
    const stk_classic::mesh::PartVector& parts = meta_data->get_parts();

    stk_classic::mesh::PartVector::const_iterator partItr;
    for ( partItr = parts.begin(); partItr != parts.end(); ++partItr ) 
    {
	const stk_classic::mesh::Part* part = *partItr;
	const stk_classic::mesh::PartVector& subsets = part->subsets();
	const CellTopologyData* ct = 
	    meta_data->get_cell_topology(*part).getCellTopologyData();

	if ( part->primary_entity_rank() == meta_data->side_rank() && 
	    ct == 0&& subsets.size() > 0 )
	{
	    TEUCHOS_TEST_FOR_EXCEPTION(
		subsets.size()!=1,std::runtime_error,
		"registerSidesets error - part \"" << 
		part->name() << "\" has more than one subset"); 

	    // grab cell topology and name of subset part
	    stk_classic::mesh::Part* ss_part = subsets[0];
	    const CellTopologyData*  ss_ct = 
		meta_data->get_cell_topology(*ss_part).getCellTopologyData();
 
	    // only add subset parts that have a topology
	    if ( ss_ct != 0 )
	    {
		sideset_parts.push_back( ss_part );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Register the nodesets in an exodus file.
void registerNodesets(
    const Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData>& meta_data, 
    stk_classic::io::MeshData& mesh_data,
    stk_classic::mesh::PartVector& nodeset_parts )
{
    nodeset_parts.clear();
    const stk_classic::mesh::PartVector& parts = meta_data->get_parts();

    stk_classic::mesh::PartVector::const_iterator partItr;
    for ( partItr = parts.begin(); partItr != parts.end(); ++partItr ) 
    {
	stk_classic::mesh::Part* part = *partItr;
	const CellTopologyData* ct = 
	    meta_data->get_cell_topology(*part).getCellTopologyData();

	// if a side part ==> this is a sideset: now storage is recursive
	// on part contains all sub parts with consistent topology
	if ( part->primary_entity_rank() == meta_data->node_rank() && ct == 0 ) 
	{
	    // only add subset parts that have no topology
	    if ( part->name() != "nodes" )
	    {
		nodeset_parts.push_back( part );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Load an exodus file into a STK database meta data.
void loadExodusFileToMetaData(
    const std::string& file_name,
    const int dimension,
    stk_classic::ParallelMachine comm,
    Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData>& meta_data,
    stk_classic::mesh::PartVector& element_block_parts,
    stk_classic::mesh::PartVector& sideset_parts,
    stk_classic::mesh::PartVector& nodeset_parts,
    stk_classic::io::MeshData& mesh_data )
{
    meta_data = Teuchos::rcp( new stk_classic::mesh::fem::FEMMetaData(dimension) );
    stk_classic::io::create_input_mesh( 
	"exodus", file_name, comm, *meta_data, mesh_data );
    registerElementBlocks( meta_data, mesh_data, element_block_parts );
    registerSidesets( meta_data, mesh_data, sideset_parts );
    registerNodesets( meta_data, mesh_data, nodeset_parts );
}

//---------------------------------------------------------------------------//
// Load an exodus file into a STK database bulk data.
void loadExodusFileToBulkData( 
    stk_classic::io::MeshData& mesh_data,
    stk_classic::ParallelMachine comm,
    Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData>& meta_data,
    Teuchos::RCP<stk_classic::mesh::BulkData>& bulk_data )
{
    meta_data->commit();
    stk_classic::io::define_input_fields( mesh_data, *meta_data );
    bulk_data = Teuchos::rcp( 
	new stk_classic::mesh::BulkData( 
	    stk_classic::mesh::fem::FEMMetaData::get_meta_data(*meta_data), comm) );
    bulk_data->modification_begin();
    stk_classic::io::populate_bulk_data( *bulk_data, mesh_data );
    bulk_data->modification_end();
}

//---------------------------------------------------------------------------//
// Write a STK database to an exodus file.
void writeExodusFile( const std::string& file_name,
		      stk_classic::ParallelMachine comm,
		      const stk_classic::mesh::fem::FEMMetaData& meta_data,
		      stk_classic::mesh::BulkData& bulk_data )
{
    stk_classic::io::MeshData mesh_data;
    stk_classic::io::create_output_mesh( file_name, comm, bulk_data, mesh_data );
    stk_classic::io::define_output_fields( mesh_data, meta_data, true );
    stk_classic::io::process_output_request( mesh_data, bulk_data, 0.0 );
}

//---------------------------------------------------------------------------//
// Evaluate the source function at a given point. 
double evaluate( const double x, const double y, 
		 const double z, const int function_type,
		 const double c, const double p )
{
    // This is the neutron flux profile of a critical finite cylinder for the
    // pin-channel problem.
    if ( 1 == function_type )
    {
	double r = std::sqrt(x*x+y*y);
	return -std::cos( 3.14159*r/3.2 )*std::cos( 3.14159*z/31.0 );
    }

    // Parabolic function for the boxes problem.
    else if ( 2 == function_type )
    {
	return 1.0 + x*x + y*y;
    }

    // Parabolic function for the cubes problem.
    else if ( 3 == function_type )
    {
	return 1.0 + x*x + y*y + z*z;
    }

    // Function for the spheres problem.
    else if ( 4 == function_type )
    {
	return 1.0 + x*y*z - x*y;
    }

    // Step function for the boxes problem.
    else if ( 5 == function_type )
    {
	return ( x >= 1.0 ) ? 2.0 : 1.0;
    }

    // Spacer grid function.
    else if ( 6 == function_type )
    {
	double r = std::sqrt((x-19)*(x-19)+(y-19)*(y-19));
	return std::cos(r/20.0) * std::cos(z/60.0);
    }

    // Second derivative function.
    else if ( 7 == function_type )
    {
	return 
	    c*( std::pow(std::abs(x/3.0),p) + 
		std::pow(std::abs(y/3.0),p) + 
		std::pow(std::abs(z/3.0),p) ) +
	    std::abs( x/5.0 ) + std::abs( y/5.0 ) + std::abs( z/5.0 );
    }

    // Default
    return 0.0;
}

//---------------------------------------------------------------------------//
// Evaluate the source function on the source nodes.
void evaluateSourceFunction( const Teuchos::ArrayView<const double>& coords,
			     const Teuchos::ArrayView<double>& data,
			     const int function_type,
			     const double c, const double p )
{
    for ( unsigned i = 0; i < data.size(); ++i )
    {
	data[i] = evaluate( coords[3*i], coords[3*i+1], coords[3*i+2],
			    function_type, c, p );
    }
}

//---------------------------------------------------------------------------//
// Source field evaluator.
class SourceEvaluator : public Topaz::FieldEvaluator
{
  public:

    typedef Intrepid::FieldContainer<double> MDArray;

    SourceEvaluator( const Teuchos::RCP<Topaz::SurfaceMesh>& surface,
		     stk_classic::mesh::Field<double>* field )
	: d_surface( surface )
	, d_mesh( surface->d_mesh )
	, d_field( field )
	, d_first_run( true )
	, d_eval_counter( 0 )
    {
	// Get the topology of the side.
	d_mesh->getEntityTopology( 
	    *d_surface->d_sides[0].d_side, d_side_topo );

	// Get the topology of the parent element.
	d_mesh->getEntityTopology( 
	    *d_surface->d_sides[0].d_parent_element, d_elem_topo );

	// Build the target basis for the element
	d_basis = Topaz::BasisFactory::create( d_elem_topo );

	// Build the side id map.
	for ( unsigned i = 0; i < d_surface->d_sides.size(); ++i )
	{
	    d_side_map.insert( 
		std::pair<stk_classic::mesh::EntityId,unsigned>(
		    d_surface->d_sides[i].d_side->identifier(),i) );
	}
    }

    ~SourceEvaluator() { /* ... */ }

    int dimension() { return 1; }
    
    // Evaluation function.
    Teuchos::Array<double> evaluate(
	const stk_classic::mesh::EntityId& entity_id,
	const Teuchos::ArrayView<const double>& evaluation_points )
    {
	// Get the face local id, parent, and side id
	unsigned local_side_id = d_side_map.find( entity_id )->second;
	stk_classic::mesh::Entity* parent = 
	    d_surface->d_sides[local_side_id].d_parent_element;
	int side_id = d_surface->d_sides[local_side_id].d_side_id;

	// Get the vertices of the parent.
	Teuchos::Array<stk_classic::mesh::Entity*> parent_nodes;
	d_mesh->getRelations( *parent, d_mesh->nodeRank(), parent_nodes );

	// Get the degrees of freedom on the parent nodes.
	Teuchos::Array<double> dofs( parent_nodes.size() );
	for ( unsigned i = 0; i < parent_nodes.size(); ++i )
	{
	    dofs[i] = *stk_classic::mesh::field_data( *d_field, *parent_nodes[i] );
	}

	// Get the canonical ids of the element nodes that are on the face.
	int num_evaluation_points = evaluation_points.size() / 3;
	int num_side_nodes = d_side_topo.getNodeCount();
	int cardinality = d_basis->getCardinality();
	Teuchos::Array<unsigned> side_nodes( num_side_nodes );
	for ( int i = 0; i < num_side_nodes; ++i )
	{
	    side_nodes[i] = d_elem_topo.getNodeMap( 2, side_id, i );
	}

	// If this is the first run, cache the basis values.
	MDArray basis_evals( cardinality, num_evaluation_points );
	if ( d_first_run )
	{
	    // Get the coordinates of the parent entity.
	    MDArray parent_coords;
	    d_mesh->getSingleEntityNodeCoordinates( *parent, parent_coords );
	    Teuchos::Array<int> cell_dims( 3, 1 );
	    cell_dims[1] = parent_coords.dimension(0);
	    cell_dims[2] = parent_coords.dimension(1);
	    MDArray cell_coords( cell_dims, parent_coords.getData()() );

	    // Map the evaluation points to the reference frame of parent.
	    Teuchos::Array<int> eval_dims( 2 );
	    eval_dims[0] = num_evaluation_points;
	    eval_dims[1] = d_mesh->getDimension();
	    MDArray physical_points( 
		eval_dims, 
		Teuchos::ArrayView<double>(
		    const_cast<double*>(evaluation_points.getRawPtr()),
		    evaluation_points.size()) );
	    MDArray reference_points( num_evaluation_points, d_mesh->getDimension() );
	    Intrepid::CellTools<double>::mapToReferenceFrame(
		reference_points, physical_points, cell_coords, d_elem_topo, 0 );

	    // Evaluate the basis in the cell.
	    d_basis->getValues( basis_evals, reference_points, Intrepid::OPERATOR_VALUE );
	    d_basis_evals.push_back(basis_evals);
	}

	// Otherwise get those that we cached.
	else
	{
	    basis_evals = d_basis_evals[d_eval_counter];
	    ++d_eval_counter;
	}

	// Extract the degree of face values and sum to give the evaluation.
	Teuchos::Array<double> evals( num_side_nodes * num_evaluation_points, 0.0 );
	for ( int i = 0; i < num_evaluation_points; ++i )
	{
	    for ( int j = 0; j < num_side_nodes; ++j )
	    {
		evals[i] += dofs[side_nodes[j]] * basis_evals(side_nodes[j],i);
	    }
	}

	// Evaluate the basis at the evaluation points.
	return evals;
    }

    Teuchos::RCP<Topaz::SurfaceMesh> d_surface;
    Teuchos::RCP<Topaz::STKMeshInterface> d_mesh;
    shards::CellTopology d_side_topo;
    shards::CellTopology d_elem_topo;
    Teuchos::RCP<Intrepid::Basis<double,MDArray> > d_basis;
    stk_classic::mesh::Field<double>* d_field;
    std::tr1::unordered_map<stk_classic::mesh::EntityId,unsigned> d_side_map;

    // This pieces are a hack to get the runtime down for some numerical
    // experiments. In general, this is bad code. Don't ever do this.
    void firstRunDone() { d_first_run = false; }
    void resetEvalCounter() { d_eval_counter = 0; }
    Teuchos::Array<MDArray> d_basis_evals;
    bool d_first_run;
    int d_eval_counter;
};

//---------------------------------------------------------------------------//
// Target basis evaluator.
class TargetEvaluator : public Topaz::FieldEvaluator
{
  public:

    typedef Intrepid::FieldContainer<double> MDArray;

    TargetEvaluator( const Teuchos::RCP<Topaz::SurfaceMesh>& surface )
	: d_surface( surface )
	, d_mesh( surface->d_mesh )
	, d_first_run( true )
	, d_eval_counter( 0 )
    { 
	// Get the topology of the side.
	d_mesh->getEntityTopology( 
	    *d_surface->d_sides[0].d_side, d_side_topo );

	// Get the topology of the parent element.
	d_mesh->getEntityTopology( 
	    *d_surface->d_sides[0].d_parent_element, d_elem_topo );

	// Build the target basis for the element
	d_basis = Topaz::BasisFactory::create( d_elem_topo );

	// Build the side id map.
	for ( unsigned i = 0; i < d_surface->d_sides.size(); ++i )
	{
	    d_side_map.insert( 
		std::pair<stk_classic::mesh::EntityId,unsigned>(
		    d_surface->d_sides[i].d_side->identifier(),i) );
	}
    }

    ~TargetEvaluator() { /* ... */ }

    int dimension() { return d_side_topo.getNodeCount(); }
    
    // Evaluation function. Return the value of the basis at the points.
    Teuchos::Array<double> evaluate(
	const stk_classic::mesh::EntityId& entity_id,
	const Teuchos::ArrayView<const double>& evaluation_points )
    {
	// Compute and cache the basis values if this is the first run.
	if ( d_first_run )
	{
	    // Get the face local id, parent, and side id
	    // Get the face local id, parent, and side id
	    unsigned local_side_id = d_side_map.find( entity_id )->second;
	    stk_classic::mesh::Entity* parent = 
		d_surface->d_sides[local_side_id].d_parent_element;
	    int side_id = d_surface->d_sides[local_side_id].d_side_id;

	    // Get the canonical ids of the element nodes that are on the face.
	    int num_side_nodes = d_side_topo.getNodeCount();
	    Teuchos::Array<unsigned> side_nodes( num_side_nodes );
	    for ( int i = 0; i < num_side_nodes; ++i )
	    {
		side_nodes[i] = d_elem_topo.getNodeMap( 2, side_id, i );
	    }

	    // Get the coordinates of the parent entity.
	    MDArray parent_coords;
	    d_mesh->getSingleEntityNodeCoordinates( *parent, parent_coords );
	    Teuchos::Array<int> cell_dims( 3, 1 );
	    cell_dims[1] = parent_coords.dimension(0);
	    cell_dims[2] = parent_coords.dimension(1);
	    MDArray cell_coords( cell_dims, parent_coords.getData()() );

	    // Map the evaluation points to the reference frame of parent.
	    int num_evaluation_points = evaluation_points.size() / 3;
	    Teuchos::Array<int> eval_dims( 2 );
	    eval_dims[0] = num_evaluation_points;
	    eval_dims[1] = d_mesh->getDimension();
	    MDArray physical_points( 
		eval_dims, 
		Teuchos::ArrayView<double>(
		    const_cast<double*>(evaluation_points.getRawPtr()),
		    evaluation_points.size()) );
	    MDArray reference_points( num_evaluation_points, d_mesh->getDimension() );
	    Intrepid::CellTools<double>::mapToReferenceFrame(
		reference_points, physical_points, cell_coords, d_elem_topo, 0 );

	    // Evaluate the basis in the cell.
	    int cardinality = d_basis->getCardinality();
	    MDArray basis_evals( cardinality, num_evaluation_points );
	    d_basis->getValues( basis_evals, reference_points, Intrepid::OPERATOR_VALUE );

	    // Extract the face values and order the basis evaluations as
	    // (POINT,CARDINALITY). Because we are on the face of the element the
	    // nodes that are not on the face will have evaluated to zero.
	    Teuchos::Array<double> evals( num_side_nodes * num_evaluation_points );
	    for ( int i = 0; i < num_evaluation_points; ++i )
	    {
		for ( int j = 0; j < num_side_nodes; ++j )
		{
		    evals[num_side_nodes*i + j] = basis_evals(side_nodes[j],i);
		}
	    }
	    d_evals.push_back(evals);
	    return evals;
	}

	// Otherwise just return the cached values.
	++d_eval_counter;
	return d_evals[d_eval_counter-1];
    }

    Teuchos::RCP<Topaz::SurfaceMesh> d_surface;
    Teuchos::RCP<Topaz::STKMeshInterface> d_mesh;
    shards::CellTopology d_side_topo;
    shards::CellTopology d_elem_topo;
    Teuchos::RCP<Intrepid::Basis<double,MDArray> > d_basis;
    std::tr1::unordered_map<stk_classic::mesh::EntityId,unsigned> d_side_map;

    // This pieces are a hack to get the runtime down for some numerical
    // experiments. In general, this is bad code. Don't ever do this.
    void firstRunDone() { d_first_run = false; }
    void resetEvalCounter() { d_eval_counter = 0; }
    void resetData() { d_evals.clear(); }
    Teuchos::Array<Teuchos::Array<double> > d_evals;
    bool d_first_run;
    int d_eval_counter;
};

//---------------------------------------------------------------------------//
// Evaluate the integral of a field over a surface.
double fieldIntegral( const Topaz::SurfaceMesh& surface,
		      const stk_classic::mesh::Field<double>& field )
{
    // Get the topology of the surface sides.
    shards::CellTopology elem_topo;
    surface.d_mesh->getEntityTopology( 
	*surface.d_sides[0].d_parent_element, elem_topo );
    shards::CellTopology side_topo;				       
    surface.d_mesh->getEntityTopology( *surface.d_sides[0].d_side, side_topo );

    // Build a basis for the parent topology.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
		 cell_basis = Topaz::BasisFactory::create( elem_topo );

    // Compute the integral in each element.
    double integral = 0.0;
    int num_sides = surface.d_sides.size();
    Intrepid::FieldContainer<double> parent_vertices(
	1, elem_topo.getNodeCount(), 3 );
    for ( int i = 0; i < num_sides; ++i )
    {
	// Extract the side node coordinates.
	surface.d_mesh->getEntityNodeCoordinates(
	    Teuchos::Array<stk_classic::mesh::Entity*>(
		1,surface.d_sides[i].d_parent_element),
	    parent_vertices );

	// Get the parent element nodes.
	Teuchos::Array<stk_classic::mesh::Entity*> elem_nodes;
	surface.d_mesh->getRelations( *surface.d_sides[i].d_parent_element, 
				      surface.d_mesh->nodeRank(), 
				      elem_nodes );
	Topaz::IntrepidSideCell side_cell( side_topo, surface.d_sides[i].d_side_id,
					   elem_topo, 5 );
	Topaz::IntrepidCell::updateState( side_cell, parent_vertices );

	// Get the side node dofs.
	Intrepid::FieldContainer<double> elem_dofs( 
	    1, parent_vertices.dimension(1) );
	for ( unsigned n = 0; n < elem_nodes.size(); ++n )
	{
	    elem_dofs(0,n) = *stk_classic::mesh::field_data( field, *elem_nodes[n] );
	}

	// Get the integration coordinates.
	Intrepid::FieldContainer<double> ip_coords(
	    side_cell.getNumIntegrationPoints(),
	    side_cell.getSpatialDimension() );
	side_cell.getReferenceIntegrationCoordinates( ip_coords );

	// Evaluate the basis at the integration coordinates.
	Intrepid::FieldContainer<double> basis_evals(
	    cell_basis->getCardinality(), side_cell.getNumIntegrationPoints() );
	cell_basis->getValues( basis_evals, ip_coords, 
			       Intrepid::OPERATOR_VALUE );

	// Collapse the basis evaluations.
	Intrepid::FieldContainer<double> dof_evals(
	    1, side_cell.getNumIntegrationPoints() );
	for ( int n = 0; n < cell_basis->getCardinality(); ++n )
	{
	    for ( int p = 0; p < side_cell.getNumIntegrationPoints(); ++p )
	    {
		dof_evals(0,p) += elem_dofs(0,n)*basis_evals(n,p);
	    }
	}
	    
	// Perform the integral.
	Intrepid::FieldContainer<double> cell_integral( 1 );
	side_cell.integrate( dof_evals, cell_integral );
	integral += cell_integral(0);
    }

    // Scale the result by the total surface area.
    return integral;
}

//---------------------------------------------------------------------------//
// Example driver.
//---------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
    // INITIALIZATION
    // --------------

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    Teuchos::RCP<Teuchos::FancyOStream>
	out = Teuchos::VerboseObjectBase::getDefaultOStream();

    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Read in command line options.
    std::string xml_input_filename;
    Teuchos::CommandLineProcessor clp(false);
    clp.setOption( "xml-in-file",
		   &xml_input_filename,
		   "The XML file to read into a parameter list" );
    clp.parse(argc,argv);

    // Build the parameter list from the xml input.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile(
	xml_input_filename, Teuchos::inoutArg(*plist) );

    // Read command-line options
    std::string source_mesh_file = plist->get<std::string>("Source Mesh File");
    std::string target_mesh_file = plist->get<std::string>("Target Mesh File");
    int degree = plist->get<int>("Integration Degree");
    int function_type = plist->get<int>("Function Type");
    int num_iters = plist->get<int>("Data Transfer Iterations");
    double c = plist->get<double>("c coefficient");
    double p = plist->get<double>("p coefficient");

    // STK I/O Init.
    Ioss::Init::Initializer io;

    // Basic problem info.
    int dimension = 3;

    // Get the raw mpi communicator (basic typedef in STK).
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    stk_classic::ParallelMachine parallel_machine = (*opaque_comm)();

    // SOURCE MESH SETUP
    //----------------

    // Load the source mesh.
    stk_classic::mesh::PartVector source_element_block_parts;
    stk_classic::mesh::PartVector source_sideset_parts;
    stk_classic::mesh::PartVector source_nodeset_parts;
    Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData> source_meta_data;
    Teuchos::RCP<stk_classic::mesh::BulkData> source_bulk_data;
    stk_classic::io::MeshData source_mesh_data;

    // Load the meta data.
    loadExodusFileToMetaData( 
	source_mesh_file, dimension, parallel_machine,
	source_meta_data, source_element_block_parts,
	source_sideset_parts, source_nodeset_parts, source_mesh_data );

    // Put a field on the source mesh nodes.
    source_meta_data->declare_field<stk_classic::mesh::Field<double> >("u");
    stk_classic::mesh::Field<double>* source_field = 
	source_meta_data->get_field<stk_classic::mesh::Field<double> >("u");
    
    // The data will be defined on the nodes of the problem.
    stk_classic::mesh::Part& source_universal = source_meta_data->universal_part();
    stk_classic::mesh::put_field( *source_field, source_meta_data->node_rank(), 
			  source_universal );

    // Add parts and fields for output.
    stk_classic::io::put_io_part_attribute( source_universal );
    stk_classic::io::set_field_role( *source_field, Ioss::Field::TRANSIENT );

    // Load the bulk data.
    loadExodusFileToBulkData( source_mesh_data, parallel_machine,
			      source_meta_data, source_bulk_data );

    // Construct a STK interface for the source mesh.
    Teuchos::RCP<Topaz::STKMeshInterface> source_mesh = Teuchos::rcp(
	new Topaz::STKMeshInterface(source_meta_data, source_bulk_data) );
    source_mesh->createAdjacentEntities();

    // Get the surface of the source mesh.
    Teuchos::RCP<Topaz::SurfaceMesh> source_surface = Teuchos::rcp(
	new Topaz::SurfaceMesh(
	    source_mesh, source_sideset_parts, source_element_block_parts) );


    // TARGET MESH SETUP
    //----------------

    // Load the target mesh.
    stk_classic::mesh::PartVector target_element_block_parts;
    stk_classic::mesh::PartVector target_sideset_parts;
    stk_classic::mesh::PartVector target_nodeset_parts;
    Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData> target_meta_data;
    Teuchos::RCP<stk_classic::mesh::BulkData> target_bulk_data;
    stk_classic::io::MeshData target_mesh_data;

    // Load the meta data.
    loadExodusFileToMetaData( 
	target_mesh_file, dimension, parallel_machine,
	target_meta_data, target_element_block_parts,
	target_sideset_parts, target_nodeset_parts, target_mesh_data );

    // Put a field on the target mesh nodes.
    target_meta_data->declare_field<stk_classic::mesh::Field<double> >("u");
    stk_classic::mesh::Field<double>* target_field = 
	target_meta_data->get_field<stk_classic::mesh::Field<double> >("u");
    
    // The data will be defined on the nodes of the problem.
    stk_classic::mesh::Part& target_universal = target_meta_data->universal_part();
    stk_classic::mesh::put_field( *target_field, target_meta_data->node_rank(), 
			  target_universal );

    // Add parts and fields for output.
    stk_classic::io::put_io_part_attribute( target_universal );
    stk_classic::io::set_field_role( *target_field, Ioss::Field::TRANSIENT );

    // Load the bulk data.
    loadExodusFileToBulkData( target_mesh_data, parallel_machine,
			      target_meta_data, target_bulk_data );

    // Construct a STK interface for the target mesh.
    Teuchos::RCP<Topaz::STKMeshInterface> target_mesh = Teuchos::rcp(
	new Topaz::STKMeshInterface(target_meta_data, target_bulk_data) );
    target_mesh->createAdjacentEntities();

    // Get the surface of the target mesh.
    Teuchos::RCP<Topaz::SurfaceMesh> target_surface = Teuchos::rcp(
	new Topaz::SurfaceMesh(
	    target_mesh, target_sideset_parts, target_element_block_parts) );


    // PERFORM THE TRANSFER
    // ------------------------------

    // Get the source mesh node coordinates on the surface.
    int num_source_nodes = source_surface->d_nodes.size();
    Teuchos::Array<double> source_surface_coords( 
	dimension*num_source_nodes );
    double* node_coords;
    std::tr1::unordered_set<stk_classic::mesh::Entity*>::const_iterator node_iterator;
    int lid = 0;
    for ( node_iterator = source_surface->d_nodes.begin();
	  node_iterator != source_surface->d_nodes.end();
	  ++node_iterator )
    {
	node_coords = source_mesh->getNodeCoordinates( **node_iterator );
	std::copy( node_coords, node_coords+dimension,
		   &source_surface_coords[lid*dimension] );
	++lid;
    }

    // Get the target mesh node coordinates on the surface.
    int num_target_nodes = target_surface->d_nodes.size();
    Teuchos::Array<double> target_surface_coords( 
	dimension*num_target_nodes );
    lid = 0;
    for ( node_iterator = target_surface->d_nodes.begin();
	  node_iterator != target_surface->d_nodes.end();
	  ++node_iterator )
    {
	node_coords = target_mesh->getNodeCoordinates( **node_iterator );
	std::copy( node_coords, node_coords+dimension,
		   &target_surface_coords[lid*dimension] );
	++lid;
    }

    // Create a source-to-target field evaluator.
    Teuchos::RCP<Topaz::FieldEvaluator> s2t_source_evaluator = 
	Teuchos::rcp( new SourceEvaluator(source_surface,source_field) );

    // Create a source-to-target target basis evaluator.
    Teuchos::RCP<Topaz::FieldEvaluator> s2t_target_evaluator = 
	Teuchos::rcp( new TargetEvaluator(target_surface) );

    // Build the source-to-target interpolation object.
    Topaz::Tolerances s2t_tolerances;
    Topaz::TransferOperator s2t_transfer_operator(
	degree, source_surface, target_surface, s2t_target_evaluator, s2t_tolerances );
    s2t_transfer_operator.setup();

    // Create a target-to-source field evaluator.
    Teuchos::RCP<Topaz::FieldEvaluator> t2s_source_evaluator = 
	Teuchos::rcp( new SourceEvaluator(target_surface,target_field) );

    // Create a target-to-source target basis evaluator.
    Teuchos::RCP<Topaz::FieldEvaluator> t2s_target_evaluator = 
	Teuchos::rcp( new TargetEvaluator(source_surface) );

    // Build the target-to-source interpolation object.
    Topaz::Tolerances t2s_tolerances;
    Topaz::TransferOperator t2s_transfer_operator(
	degree, target_surface, source_surface, t2s_target_evaluator, t2s_tolerances );
    t2s_transfer_operator.setup();

    // Evaluate the function at the source mesh nodes.
    Teuchos::ArrayRCP<double> source_data( num_source_nodes );
    evaluateSourceFunction( source_surface_coords(), source_data(), function_type, c, p );

    // Allocate space for the target mesh data.
    Teuchos::ArrayRCP<double> target_data( num_target_nodes );

    // Compute source and target reference solutions.
    Teuchos::Array<double> source_gold( num_source_nodes );
    evaluateSourceFunction( source_surface_coords(), source_gold(), function_type, c, p );
    double source_gold_l2 = 0.0;
    for ( int i = 0; i < num_source_nodes; ++i )
    {
	source_gold_l2 += source_gold[i]*source_gold[i];
    }
    source_gold_l2 = std::sqrt( source_gold_l2 );
    Teuchos::Array<double> target_gold( num_target_nodes );
    evaluateSourceFunction( target_surface_coords(), target_gold(), function_type, c, p );
    double target_gold_l2 = 0.0;
    for ( int i = 0; i < num_target_nodes; ++i )
    {
	target_gold_l2 += target_gold[i]*target_gold[i];
    }
    target_gold_l2 = std::sqrt( target_gold_l2 );

    // Apply the initial data to the source mesh.
    lid = 0;
    for ( node_iterator = source_surface->d_nodes.begin();
	  node_iterator != source_surface->d_nodes.end();
	  ++node_iterator )
    {
	*stk_classic::mesh::field_data( *source_field, **node_iterator ) = 
	    source_data[lid];
	++lid;
    }

    // Integrate the field over the source surface to get the initial value.
    double gold_integral = fieldIntegral( *source_surface, *source_field );

    // Perform data transfer iterations. Check accuracy and convergence at
    // each iteration.
    Teuchos::rcp_dynamic_cast<TargetEvaluator>(t2s_target_evaluator)->resetData();
    Teuchos::rcp_dynamic_cast<TargetEvaluator>(s2t_target_evaluator)->resetData();
    std::cout << std::endl;
    std::cout << "Iteration   ||s||          ||t||            S              T           " << std::endl;
    std::cout << "=======================================================================" << std::endl;
    for ( int i = 0; i < num_iters; ++i )
    {
	// Transfer the data from the source to the target.
	Teuchos::rcp_dynamic_cast<SourceEvaluator>(s2t_source_evaluator)->resetEvalCounter();
	Teuchos::rcp_dynamic_cast<TargetEvaluator>(s2t_target_evaluator)->resetEvalCounter();
	s2t_transfer_operator.apply( s2t_source_evaluator, target_data );
	if ( i == 0 )
	{
	    Teuchos::rcp_dynamic_cast<SourceEvaluator>(s2t_source_evaluator)->firstRunDone();
	    Teuchos::rcp_dynamic_cast<TargetEvaluator>(s2t_target_evaluator)->firstRunDone();
	}

	// Compute the target accuracy error and apply the data to the mesh.
	double local_error = 0.0;
	double target_l2_accuracy = 0.0;
	lid = 0;
	for ( node_iterator = target_surface->d_nodes.begin();
	      node_iterator != target_surface->d_nodes.end();
	      ++node_iterator )
	{
	    local_error = (target_data[lid]-target_gold[lid]);
	    target_l2_accuracy += local_error*local_error;

	    *stk_classic::mesh::field_data( *target_field, **node_iterator ) = 
		target_data[lid];

	    ++lid;
	}

	// Integrate the field over the target surface.
	double target_integral = fieldIntegral( *target_surface, *target_field );

	// Compute the target conservation error.
	double target_conservation_error =
	    std::abs((target_integral - gold_integral) / gold_integral);

	// Transfer the data from the source to the target.
	Teuchos::rcp_dynamic_cast<SourceEvaluator>(t2s_source_evaluator)->resetEvalCounter();
	Teuchos::rcp_dynamic_cast<TargetEvaluator>(t2s_target_evaluator)->resetEvalCounter();
	t2s_transfer_operator.apply( t2s_source_evaluator, source_data );
	if ( i == 0 )
	{
	    Teuchos::rcp_dynamic_cast<SourceEvaluator>(t2s_source_evaluator)->firstRunDone();
	    Teuchos::rcp_dynamic_cast<TargetEvaluator>(t2s_target_evaluator)->firstRunDone();
	}

	// Compute the source accuracy error and apply the data to the mesh.
	double source_l2_accuracy = 0.0;
	lid = 0;
	for ( node_iterator = source_surface->d_nodes.begin();
	      node_iterator != source_surface->d_nodes.end();
	      ++node_iterator )
	{
	    local_error = (source_data[lid]-source_gold[lid]);
	    source_l2_accuracy += local_error*local_error;

	    *stk_classic::mesh::field_data( *source_field, **node_iterator ) = 
		source_data[lid];
	    ++lid;
	}

	// Integrate the field over the source surface.
	double source_integral = fieldIntegral( *source_surface, *source_field );

	// Compute the source conservation error.
	double source_conservation_error = 
	    std::abs((source_integral - gold_integral) / gold_integral);

	// Output the results.
	std::cout << i << "          " 
		  << std::sqrt(source_l2_accuracy) / source_gold_l2 << "   "
		  << std::sqrt(target_l2_accuracy) / target_gold_l2 << "   "
		  << source_conservation_error << "    "
		  << target_conservation_error << std::endl;
    }
    std::cout << std::endl;

    // WRITE THE MESHES TO FILE
    // -------------------------
    std::string source_out_file = "out_" + source_mesh_file;
    writeExodusFile( source_out_file, parallel_machine, 
		     *source_meta_data, *source_bulk_data );

    std::string target_out_file = "out_" + target_mesh_file;
    writeExodusFile( target_out_file, parallel_machine, 
		     *target_meta_data, *target_bulk_data );
}

//---------------------------------------------------------------------------//
// end tstSTK_Mesh.cpp
//---------------------------------------------------------------------------//
