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

#include "DTK_FieldTraits.hpp"
#include <DTK_MeshTraits.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class DataSource
 * \brief DataSource objects provide mesh and methods to evaluate fields
 * associated with that mesh.
 */
//===========================================================================//
template<class Mesh, class DataField>
class DataSource
{
  public:

    // Forward declare the field evaluator.
    class FieldEvaluator;

    //@{
    //! Typedefs.
    typedef Mesh                             mesh_type;
    typedef MeshTraits<Mesh>                 MT;
    typedef MT::handle_type                  handle_type;
    typedef DataField                        data_field_type;
    typedef FieldTraits<DataField>           FT;
    typedef FieldEvaluator<handle_type>      FieldEvaluatorType;
    typedef Teuchos::RCP<FieldEvaluatorType> RCP_FieldEvaluator;
    typedef Teuchos::Comm<int>               CommType;
    typedef Teuchos::RCP<CommType>           RCP_Comm;
    //@}


    /*!
     * \brief Mix-in interface for field evaluation kernels. 
     *
     * This is inside the class definition to ensure type consistency.
     */
    class FieldEvaluator
    {
      public:

	//@{
	//! Typedefs.
	typedef typename MT::handle_type     handle_type;
	//@}

	//! Constructor.
	FieldEvaluator()
	{ /* ... */ }

	//! Destructor.
	virtual ~FieldEvaluator()
	{ /* ... */ }

	/*!
	 * \brief Evaluate the field in the given elements at the given
	 * coordinates and return the evaluations in a DataField.
	 * \param elements A vector of element handles in which to evaluate
	 * the field.
	 * \param coords A vector of interleaved coordinates 
	 * ( x0, y0, z0, ..., xN, yN, zN ) at which to evaluate the
	 * field. Coordinates ( x_n, y_n, z_n ) should be evaluated in the nth
	 * element in the elements vector.
	 * \return Return a DataField containing the evaluated field
	 * values. This returned field is expected to be of the same length as
	 * the elements input vector. For those coordinates that can't be
	 * evaluated in the given element, return 0 in their position.
	 */
	virtual DataField evaluate( const std::vector<handle_type>& elements,
				    const std::vector<double>& coords ) = 0;
    };


  public:

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
