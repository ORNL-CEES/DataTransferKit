//---------------------------------------------------------------------------//
/*!
 * \file DTK_DataTraits.hpp
 * \author Stuart R. Slattery
 * \brief Traits declaration for field types.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATAFIELDTRAITS_HPP
#define DTK_DATAFIELDTRAITS_HPP

#include <iterator>

namespace DataTransferKit
{

/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename T>
struct UndefinedFieldTraits
{
    static inline T notDefined() 
    { return T::this_type_is_missing_a_specialization(); }
};

/*!
 * \brief Field traits definitions.
 * 
 * These traits correlate to the basic concept of a field within DTK. A field
 * can contain anything in an array, but it must store its objects in
 * contiguous memory. 
 */
template<typename T>
struct FieldTraits
{
    //! Typedef for value type. The field type must implement
    //! Teuchos::ScalarTraits.
    typedef typename T::value_type value_type;
    
    //{@
    //! Typedef for field value random access iterators.
    typedef typename 
    std::iterator<std::random_access_iterator_tag,value_type> iterator;

    typedef typename 
    std::const_iterator<std::random_access_iterator_tag,value_type> const_iterator;
    //@}

    //! Returns the number of elements in the field.
    static inline std::size_t size( const T &field )
    { return UndefinedFieldTraits<T>::notDefined(); }

    //@{
    //! Returns the iterator to the front of the field.
    static inline iterator begin( T &field )
    { return UndefinedFieldTraits<T>::notDefined(); }

    static inline const_iterator begin( const T &field )
    { return UndefinedFieldTraits<T>::notDefined(); }
    //@}

    //@{
    //! Returns the iterator to the end of the field.
    static inline iterator end( T &field )
    { return UndefinedFieldTraits<T>::notDefined(); }

    static inline const_iterator end( const T &field )
    { return UndefinedFieldTraits<T>::notDefined(); }
    //@}

    //! Returns if the field is empty.
    static inline bool empty( const T &field )
    { return UndefinedFieldTraits<T>::notDefined(); }
};

} // end namespace DataTransferKit

#endif // end DTK_DATAFIELDTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_DataFieldTraits.hpp
//---------------------------------------------------------------------------//

