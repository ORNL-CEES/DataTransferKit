//---------------------------------------------------------------------------//
/*!
 * \file   DTK_DataSource.hpp
 * \author Stuart R. Slattery
 * \brief  DataSource declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATASOURCE_HPP
#define DTK_DATASOURCE_HPP

#include <string>
#include <vector>
#include <map>

#include <mpi.h>

#include "DTK_FieldEvaluator.hpp"
#include "DTK_FieldTraits.hpp"
#include <DTK_MeshTraits.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class DataSource
 * \brief DataSource objects provide mesh and kernels to evaluate fields
 * associated with that mesh.
 */
//---------------------------------------------------------------------------//
template<class Mesh, class DataField>
class DataSource
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                mesh_type;
    typedef MeshTraits<Mesh>                    MT;
    typedef typename MT::global_ordinal_type    global_ordinal_type;
    typedef DataField                           data_field_type;
    typedef FieldTraits<DataField>              FT;
    typedef typename FT::value_type             data_type;
    typedef FieldEvaluator<Mesh,DataField>      FieldEvaluatorType;
    typedef Teuchos::RCP<FieldEvaluatorType>    RCP_FieldEvaluator;
    typedef Teuchos::Comm<int>                  CommType;
    typedef Teuchos::RCP<CommType>              RCP_Comm;
    //@}

    // Constructor.
    DataSource( const Mesh& mesh, const MPI_Comm& mpi_comm );

    // Destructor.
    ~DataSource();

    // Register a source field evaluator.
    void registerFieldEvaluator( const std::string& field_name,
				 const RCP_FieldEvaluator& field_evaluator );

    //! Get the id for a field given its name.
    std::size_t getFieldId( const std::string& name ) const
    { return d_name_map.find( name )->second; }

    //! Get the evaluation kernel for a field given its id.
    const RCP_FieldEvaluator& getFieldEvaluator( const std::size_t id ) const
    { return d_eval_map.find( id )->second; }
    
    //! Get the mesh.
    const Mesh& getMesh() const
    { return d_mesh; }

    //! Get the communicator.
    const RCP_Comm& getComm() const
    { return d_comm; }

  private:

    // Mesh.
    Mesh d_mesh;
   
    // Communicator.
    RCP_Comm d_comm;

    // Name to id map.
    std::map<std::string,std::size_t> d_name_map;

    // Id to evaluation map.
    std::map<std::size_t,RCP_FieldEvaluator> d_eval_map;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_DataSource_def.hpp"

#endif // DTK_DATASOURCE_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataSource.hpp
//---------------------------------------------------------------------------//
