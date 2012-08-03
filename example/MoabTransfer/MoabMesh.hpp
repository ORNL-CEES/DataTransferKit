//---------------------------------------------------------------------------//
/*!
 * \file MoabMesh.hpp
 * \author Stuart R. Slattery
 * \brief Moab mesh declaration for example.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABMESH_EX_HPP
#define DTK_MOABMESH_EX_HPP

#include <string>
#include <vector>

#include "ArrayField.hpp"

#include <DTK_MeshTraits.hpp>

#include <MBInterface.hpp>
#include <MBRange.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

//---------------------------------------------------------------------------//
class MoabMesh
{
  public:

    typedef moab::EntityHandle                        global_ordinal_type;
    typedef Teuchos::Comm<int>                        CommType;
    typedef Teuchos::RCP<const CommType>              RCP_Comm;

    MoabMesh()
    { /* ... */ }

    MoabMesh( const RCP_Comm& comm, 
	      const std::string& filename, 
	      const moab::EntityType& block_topology,
	      const int partitioning_type );

    ~MoabMesh()
    { /* ... */ }

    const Teuchos::RCP<moab::Interface>& getMoab() const
    { return d_moab; }
    
    void tag( const ArrayField& data );

    void write( const std::string& filename );

    std::vector<global_ordinal_type>::const_iterator nodesBegin() const
    { return d_vertices.begin(); }

    std::vector<global_ordinal_type>::const_iterator nodesEnd() const
    { return d_vertices.end(); }

    std::vector<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    std::vector<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    const std::vector<double>& getCoords() const
    { return d_coords; }

    std::size_t nodeDim() const
    { return d_node_dim; }

    std::size_t blockTopology() const;

    std::size_t nodesPerElement() const
    { return d_nodes_per_element; }

    std::vector<global_ordinal_type>::const_iterator elementsBegin() const
    { return d_elements.begin(); }

    std::vector<global_ordinal_type>::const_iterator elementsEnd() const
    { return d_elements.end(); }

    std::vector<global_ordinal_type>::const_iterator connectivityBegin() const
    { return d_connectivity.begin(); }

    std::vector<global_ordinal_type>::const_iterator connectivityEnd() const
    { return d_connectivity.end(); }
    
    std::vector<std::size_t>::const_iterator permutationBegin() const
    { return d_permutation_list.begin(); }

    std::vector<std::size_t>::const_iterator permutationEnd() const
    { return d_permutation_list.end(); }

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Moab interface.
    Teuchos::RCP<moab::Interface> d_moab;

    // Block topology.
    moab::EntityType d_topology;

    // Moab mesh vertices.
    std::vector<global_ordinal_type> d_vertices;

    // Coordinate blocks.
    std::vector<double> d_coords;

    // Node dimension.
    std::size_t d_node_dim;

    // Nodes per element.
    std::size_t d_nodes_per_element;

    // Moab mesh elements.
    std::vector<global_ordinal_type> d_elements;

    // Element connectivity.
    std::vector<global_ordinal_type> d_connectivity;

    // Permutation list.
    std::vector<std::size_t> d_permutation_list;
};

//---------------------------------------------------------------------------//
// Mesh traits definition.
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

template<>
class MeshTraits<MoabMesh>
{
  public:

    typedef MoabMesh::global_ordinal_type global_ordinal_type;
    typedef std::vector<global_ordinal_type>::const_iterator 
    const_node_iterator;
    typedef std::vector<double>::const_iterator 
    const_coordinate_iterator;
    typedef std::vector<global_ordinal_type>::const_iterator 
    const_element_iterator;
    typedef std::vector<global_ordinal_type>::const_iterator 
    const_connectivity_iterator;
    typedef std::vector<std::size_t>::const_iterator 
    const_permutation_iterator;

    static inline std::size_t nodeDim( const MoabMesh& mesh )
    { return mesh.nodeDim(); }

    static inline const_node_iterator nodesBegin( const MoabMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MoabMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MoabMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MoabMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementTopology( const MoabMesh& mesh )
    { return mesh.blockTopology(); }

    static inline std::size_t nodesPerElement( const MoabMesh& mesh )
    { return mesh.nodesPerElement(); }


    static inline const_element_iterator elementsBegin( const MoabMesh& mesh )
    { return mesh.elementsBegin(); }

    static inline const_element_iterator elementsEnd( const MoabMesh& mesh )
    { return mesh.elementsEnd(); }

    static inline const_connectivity_iterator 
    connectivityBegin( const MoabMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator 
    connectivityEnd( const MoabMesh& mesh )
    { return mesh.connectivityEnd(); }

    static inline const_permutation_iterator 
    permutationBegin( const MoabMesh& mesh )
    { return mesh.permutationBegin(); }

    static inline const_permutation_iterator 
    permutationEnd( const MoabMesh& mesh )
    { return mesh.permutationEnd(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABMESH_EX_HPP

//---------------------------------------------------------------------------//
// end MoabMesh.hpp
//---------------------------------------------------------------------------//

