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

#include <mpi.h>

#include "DTK_MeshTraits.hpp"
#include "DTK_FieldTraits.hpp"

#include <Teuchos_Describable.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class DataSource
 * \brief Protocol declaration for applications acting as a data source in
 * multiphysics coupling. DTK will neither allocate or deallocate memory in
 * the source application through this interface.
 *
 * Implementing this interface will then provide a single mesh and the fields
 * associated with that mesh. Multiple meshes will need multiple data source
 * interfaces specified.
 */
/*! 
 * \example core/test/tstInterfaces.cc
 *
 * Test of DataSource.
 */
//===========================================================================//
template<class Mesh, class DataField>
class DataSource : Teuchos::Describable
{
  public:

    //@{
    //! Mesh typedefs.
    typedef Mesh                                        MeshType;
    typedef typename MeshTraits<Mesh>::handle_type      handle_type;
    typedef std::vector<handle_type>                    HandleVector;
    typedef std::vector<double>                         CoordinateVector;
    //@}

    //@{
    //! Data typedefs.
    typedef DataField                                   data_field_type;
    typedef typename FieldTraits<DataField>::value_type data_type;
    //@}

    /*!
     * \brief Constructor.
     */
    DataSource()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~DataSource()
    { /* ... */ }

    /*!
     * \brief Get the communicator object for the physics implementing this
     * interface.
     * \param source_comm The MPI communicator for the source application.
     */
    virtual const MPI_Comm& getSourceComm() = 0;

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported.
     * \param field_name The name of the field for which support is being
     * checked.
     */
    virtual bool isFieldSupported( const std::string &field_name ) = 0;

    /*! 
     * \brief Provide the local source mesh.
     * \return A const reference to a type that implements mesh traits.
     */
    virtual const Mesh& getSourceMesh() = 0;

    /*! 
     * \brief Provide a copy of the local source data evaluated in the
     * given source mesh elements at the given nodes.
     * \param field_name The name of the field to provide data from.
     * \param element_handles The local element handles in which to evaluate
     * the field. See typedef for type construction.
     * \param coordiantes The coordinates at which to evaluate the field. This
     * vector will be the 3 times the length of the element handle vector and
     * be returned in an interleaved order ( x0, y0, z0, ..., xN, yN, zN
     * ). See typedef for type construction.
     * \return The evaluated field data at the nodes. The DataField type
     * is expected to implement FieldTraits. DataField::value_type is
     * expected to implement Teuchos::ScalarTraits. 
     */
    virtual const DataField
    evaluateFieldOnTargetNodes( const std::string& field_name,
				const HandleVector& element_handles,
				const CoordinateVector& coordinates ) = 0;
};

} // end namespace DataTransferKit

#endif // DTK_DATASOURCE_HPP

//---------------------------------------------------------------------------//
// end of DTK_DataSource.hpp
//---------------------------------------------------------------------------//
