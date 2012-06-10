//---------------------------------------------------------------------------//
/*!
 * \file   DTK_DataSource.hpp
 * \author Stuart R. Slattery
 * \brief  Interface declaration for data source applications.
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
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

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

    //@{
    //! Typedefs.
    typedef Mesh                             mesh_type;
    typedef MeshTraits<Mesh>                 MT;
    typedef DataField                        data_field_type;
    typedef FieldTraits<DataField>           FT;
    typedef Teuchos::Comm<int>               CommType;
    typedef Teuchos::RCP<CommType>           RCP_Comm;
    //@}


    /*!
     * \brief Mix-in interface for field evaluation kernels. This is inside
     * the class definition to ensure type consistency.
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
	 * \param coords a vector of interleaved coordinates 
	 * ( x0, y0, z0, ..., xN, yN, zN ) at which to evaluate the
	 * field. Coordinates ( x_n, y_n, z_n ) should be evaluated in the nth
	 * element in the elements vector.
	 * \return Return a DataField containing the evaluated field
	 * values. This returned field is expected to be of the same length as
	 * the elements input vector. For those coordinates that can't be
	 * evaluated in the given element, return 0.
	 */
	virtual DataField evaluate( const std::vector<handle_type> elements,
				    const std::vector<double>& coords ) = 0;
    };


  public:

    // Constructor.
    DataSource( const MPI_Comm& mpi_comm, const Mesh& mesh );

    // Destructor.
    ~DataSource();

    //! Get the id for a field given its name.
    std::size_t getFieldId( const std::string& name ) const
    { return d_name_map.find( name )->second; }

    //! Get the evaluation kernel for a field given its id.
    FieldEvaluator& getFieldEvaluator( const std::size_t id ) const
    { return d_field_map( id )->second; }
    
    //! Get the communicator.
    const RCP_Comm& getComm() const
    { return d_comm; }

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Name to id map.
    std::map<std::string,std::size_t> d_name_map;

    // Id to evaluation map.
    std::map<std::size_t,FieldEvaluator> d_field_map;

};

} // end namespace DataTransferKit

#endif // DTK_DATASOURCE_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataSource.hpp
//---------------------------------------------------------------------------//
