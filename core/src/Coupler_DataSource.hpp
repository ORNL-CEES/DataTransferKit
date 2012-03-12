//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Coupler_DataSource.hpp
 * \author Stuart Slattery
 * \date   Thu Nov 17 07:53:43 2011
 * \brief  Interface declaration for data source applications.
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_DATASOURCE_HPP
#define COUPLER_DATASOURCE_HPP

#include <string>

#include "Coupler_Point.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Describable.hpp>

namespace Coupler
{

//===========================================================================//
/*!
 * \class DataSource
 * \brief Protocol declaration for applications acting as a data source in
 * multiphysics coupling. 
 *
 * This interface is templated on the type of field data being
 * transferred, the handle type for mesh entities, and the coordinate
 * type. Ordinal type for communication is int. 
 */
/*! 
 * \example core/test/tstInterfaces.cc
 *
 * Test of DataSource.
 */
//===========================================================================//
template<class DataType, class HandleType, class CoordinateType, int DIM>
class DataSource : Teuchos::Describable
{
  public:

    //@{
    //! Useful typedefs.
    typedef int                                        OrdinalType;
    typedef Point<DIM,HandleType,CoordinateType>       PointType;
    typedef Teuchos::Comm<OrdinalType>                 Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>         RCP_Communicator;
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
     * \return The communicator for this physics.
     */
    virtual RCP_Communicator get_source_comm() = 0;

    /*!
     * \brief Check whether or not a field is supported. Return false if this
     * field is not supported. 
     * \param field_name The name of the field for which support is being
     * checked.
     */
    virtual bool is_field_supported( const std::string &field_name ) = 0;

    /*! 
     * \brief Given (x,y,z) coordinates and an associated globally unique
     * handle encapsulated in a point object, return true if the point is in
     * the local domain, false if not. 
     * \param point Point to query the local domain with.
     */
    virtual bool is_local_point( const PointType &point ) = 0;

    /*!
     * \brief Given a list of point objects, return a list of true/false
     * values that designate if each point is in the local domain.
     * \param points View of points to query the local domain with.
     * \return Array of booleans defined to be ordered implicitly as the input
     * array. Return true if a point is in the local domain, false if not.
     */
    virtual const Teuchos::ArrayRCP<bool> 
    are_local_points( const Teuchos::ArrayView<PointType> points );

    /*! 
     * \brief Provide a const view of the local source data at the target
     * points found by is_local_point.
     * \param field_name The name of the field to provide data from.
     * \return A const view of data to be sent. There are two requirements for
     * this view: 1) it is of size equal to the number of points in the local
     * domain, 2) the data is in the same order as the points found by
     * is_local_point. This view is not required to persist as it is
     * immediately copied.
     */
    virtual const Teuchos::ArrayView<DataType> 
    get_source_data( const std::string &field_name ) = 0;

    /*!
     * \brief Given a field, get a global data element to be sent to the
     * target. 
     * \param field_name The name of the field to get data from.
     * \return The global data element.
     */
    virtual DataType 
    get_global_source_data( const std::string &field_name ) = 0;
};

} // end namespace Coupler

#include "Coupler_DataSource_Def.hpp"

#endif // COUPLER_DATASOURCE_HPP

//---------------------------------------------------------------------------//
//              end of Coupler_DataSource.hpp
//---------------------------------------------------------------------------//
