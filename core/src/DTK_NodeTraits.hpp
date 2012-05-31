//---------------------------------------------------------------------------//
/*!
 * \file DTK_NodeTraits.hpp
 * \author Stuart R. Slattery
 * \brief Declaration of node traits.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_NODESETTRAITS_HPP
#define DTK_NODESETTRAITS_HPP

#include <iterator>

namespace DataTransferKit
{

/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename T>
struct UndefinedNodeTraits
{
    static inline T notDefined() 
    { return T::this_type_is_missing_a_specialization(); }
};

/*!
 * \brief Node traits definitions.
 *
 * These traits correlate to the basic concept of a mesh node within DTK. A
 * node will consist of a globally unique handle of a type that implements
 * Teuchos::OrdinalTraits and a set of coordinates of a type that implements
 * Teuchos::ScalarTraits. Coordinates must be stored in contiguous memory.
 */
template<typename T>
struct NodeTraits
{
    //! Typedef for handle type. This type must implement
    //! Teuchos::OrdinalTraits.
    typedef typename T::handle_type handle_type;

    //! Typedef for coordinate type. This type must implement
    //! Teuchos::ScalarTraits. 
    typedef typename T::coordinate_type coordinate_type;

    //! Typedef for coordinate iterator.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, coordinate_type>  
    const_coordinate_iterator;

    //! Returns the spatial dimension of the node.
    static inline std::size_t dim()
    { return UndefinedNodeTraits<T>::notDefined(); }

    //! Returns the handle of the node.
    static inline handle_type handle( const T &node )
    { return UndefinedNodeTraits<T>::notDefined(); }

    //! Returns the iterator to the front of the coordinate array of the
    //! node. 
    static inline const_coordinate_iterator coordsBegin( const T &node )
    { return UndefinedNodeTraits<T>::notDefined(); }

    //! Returns the iterator to the end of the coordinate array of the node.
    static inline const_coordinate_iterator coordsEnd( const T &node )
    { return UndefinedNodeTraits<T>::notDefined(); }
};

} // end namespace DataTransferKit

#endif // end DTK_NODESETTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_NodeTraits.hpp
//---------------------------------------------------------------------------//

