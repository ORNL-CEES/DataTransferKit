//---------------------------------------------------------------------------//
/*!
 * \file DTK_ElementTraits.hpp
 * \author Stuart R. Slattery
 * \brief Declaration of element traits.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ELEMENTTRAITS_HPP
#define DTK_ELEMENTTRAITS_HPP

namespace DataTransferKit
{

/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename T>
struct UndefinedElementTraits
{
    static inline T notDefined() 
    { return T::this_type_is_missing_a_specialization(); }
};

/*!
 * \brief Element traits definitions.
 *
 * These traits correlate to the basic concept of a mesh element within DTK. A
 * element will consist of a globally unique handle of a type that implements
 * Teuchos::OrdinalTraits and a set of node handles of a type that implements
 * Teuchos::OrdinalTraits that coorelates to the connectivity for this
 * element. Connectivity node handles must be stored in contiguous memory.
 */
template<typename T>
struct ElementTraits
{
    //! Typedef for handle type. This type must implement
    //! Teuchos::OrdinalTraits. 
    typedef typename T::handle_type handle_type;

    //! Typdef for element connectivity random access iterator.
    typedef typename 
    std::iterator<std::random_access_iterator_tag,handle_type> 
    const_connectivity_iterator;

    //! Returns the type of the element (DTK_ElementType enum).
    static inline std::size_t type()
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Returns the topology of the element (DTK_ElementTopolpogy enum).
    static inline std::size_t topology()
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Returns the number of nodes that construct the element.
    static inline std::size_t numNodes()
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Returns the handle of the element.
    static inline handle_type handle( const T &element )
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Returns the iterator to the front of the connectivity array of the
    //! element. 
    static inline const_connectivity_iterator 
    connectivityBegin( const T &element )
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Return the iterator to the end of the connectivity array of the
    //! element. 
    static inline const_connectivity_iterator 
    connectivityEnd( const T &element )
    { return UndefinedElementTraits<T>::notDefined(); }
};

} // end namespace DataTransferKit

#endif // end DTK_ELEMENTTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_ElementTraits.hpp
//---------------------------------------------------------------------------//
