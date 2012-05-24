//---------------------------------------------------------------------------//
/*!
 * \file DTK_ElementTraits.hpp
 * \author Stuart R. Slattery
 * \brief Declaration of element traits.
 */

#ifndef DTK_NODESETTRAITS_HPP
#define DTK_NODESETTRAITS_HPP

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
    //! Required typedef for handle type. This type must implement
    //! Teuchos::OrdinalTraits. 
    typedef T::handle_type handle_type;

    //! Returns the topology of the element (DTK_ElementTopolpogy enum).
    static inline std::size_t topology()
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Returns the handle.
    static inline handle_type handle()
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Returns the iterator to the front of the connectivity array. 
    static inline handle_type* connectivityBegin()
    { return UndefinedElementTraits<T>::notDefined(); }

    //! Return the iterator to the end of the connectivity array.
    static inline handle_type* connectivityEnd()
    { return UndefinedElementTraits<T>::notDefined(); }
};

} // end namespace DataTransferKit

#endif // end DTK_ELEMENTSETTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_ElementTraits.hpp
//---------------------------------------------------------------------------//
