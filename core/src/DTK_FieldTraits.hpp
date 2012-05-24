//---------------------------------------------------------------------------//
/*!
 * \file DTK_DataTraits.hpp
 * \author Stuart R. Slattery
 * \brief Traits declaration for field types.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATAFIELDTRAITS_HPP
#define DTK_DATAFIELDTRAITS_HPP

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
    typedef T::value_type value_type;

    //! Returns the size of the field.
    static inline std::size_t size()
    { return UndefinedFieldTraits<T>::notDefined(); }

    //! Returns the iterator to the front of the field.
    static inline value_type* begin()
    { return UndefinedFieldTraits<T>::notDefined(); }

    //! Returns the iterator to the end of the field.
    static inline value_type* end()
    { return UndefinedFieldTraits<T>::notDefined(); }
}

} // end namespace DataTransferKit

#endif // end DTK_DATAFIELDTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_DataFieldTraits.hpp
//---------------------------------------------------------------------------//

