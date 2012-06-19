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
template<typename UndefinedFieldType>
struct UndefinedFieldTraits
{
    static inline UndefinedFieldType notDefined() 
    { return UndefinedFieldType::this_type_is_missing_a_specialization(); }
};

/*!
 * \brief Field traits definitions.
 * 
 * These traits correlate to the basic concept of a field within DTK. A field
 * can contain anything in an array, but it must store its objects in
 * contiguous memory that is blocked by dimension. 
 */
template<typename FieldType>
class FieldTraits
{
  public:

    //! Typedef for field type.
    typedef FieldType field_type;

    //! Typedef for value type. The field type must implement
    //! Teuchos::ScalarTraits.
    typedef typename FieldType::value_type value_type;

    //! Typedef for size type.
    typedef typename FieldType::size_type size_type;
    
    //{@
    //! Typedef for field value random access iterators.
    typedef typename 
    std::iterator<std::random_access_iterator_tag,value_type> 
    iterator;

    typedef typename 
    std::iterator<std::random_access_iterator_tag,const value_type> 
    const_iterator;
    //@}

    /*! 
     * \brief Returns the dimensionality of the field ( i.e. 1 for a scalar, 3
     * for a 3 vector, 9 for a 3x3 tensor, etc. ).
     */
    static inline std::size_t dim( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    /*! 
     * \brief Returns the number of elements in the field.
     */
    static inline size_type size( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    /*! 
     * \brief Returns if the field is empty.
     */
    static inline bool empty( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    //@{
    /*! Returns the iterator to the beginning of the field. The data is
     * required to be blocked by dimensions ( d0, d1, d2, ... , dM ) as
     * ( d0_0, d0_1, ... , d0_N, d1_0, d1_1, ... , d1_N, ... , dM_N )
     */
    static inline iterator begin( FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    static inline const_iterator begin( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }
    //@}

    //@{
    /*! 
     * \brief Returns the iterator to the end of the field. The data is
     * required to be blocked by dimensions ( d0, d1, d2, ... , dM ) as
     * ( d0_0, d0_1, ... , d0_N, d1_0, d1_1, ... , d1_N, ... , dM_N )
     */
    static inline iterator end( FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }

    static inline const_iterator end( const FieldType& field )
    { return UndefinedFieldTraits<FieldType>::notDefined(); }
    //@}
};

} // end namespace DataTransferKit

#endif // end DTK_DATAFIELDTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_DataFieldTraits.hpp
//---------------------------------------------------------------------------//

