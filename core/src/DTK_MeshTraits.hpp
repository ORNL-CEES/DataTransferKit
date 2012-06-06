//---------------------------------------------------------------------------//
/*!
 * \brief DTK_MeshTraits.hpp
 * \author Stuart R. Slattery
 * \brief Declaration of mesh traits.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTRAITS_HPP
#define DTK_MESHTRAITS_HPP

namespace DataTransferKit
{

/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename T>
struct UndefinedMeshTraits
{
    static inline T notDefined() 
    { return T::this_type_is_missing_a_specialization(); }
};

/*!
 * \brief Mesh traits definitions.
 *
 * These traits correlate to the basic concept of a mesh within DTK. A
 * mesh will consist of a globally unique list of node handles of a type that
 * implements Teuchos::OrdinalTraits ( already implemented for common ordinal
 * types ) and a set of globally unique element handles of the same
 * type. Nodes are described by a coordinate field with coordinates of type
 * double. Elements are described by a list of node handles that designate
 * their connectivity.
 */
template<typename T>
struct MeshTraits
{
    //@{
    //! Typedefs.
    //! Typedef for handle type. This type must implement
    //! Teuchos::OrdinalTraits.
    typedef typename T::handle_type handle_type;

    //! Typedef for random access const iterators to handle type values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, handle_type>  
    const_handle_iterator;

    //! Typedef for random access const iterators to coordinate values. This
    //! is enforcing a coordinate type of double.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, double>  
    const_coordinate_iterator;
    //@}

    //@{
    //! Mesh node concepts.
    /*!
     * \brief Return the const iterator to the beginning of the node handle
     * block in this mesh. 
    */
    static inline const_handle_iterator nodesBegin( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*!
     * \brief Return the const iterator to the end of the node handle block in
     * this mesh.
    */ 
    static inline const_handle_iterator nodesEnd( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*!
     * \brief Return true if the coordinate block in this mesh is interleaved
     * ( x0, y0, z0, ..., xN, yN, zN ) and false if blocked 
     * ( x0, x1, ... , xN, y0, y1, ..., yN, z0, z1, ... zN ).
     */    
    static inline bool interleavedCoordinates( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*!
     * \brief Return the const iterator to the beginning of the node
     * coordinate block in this mesh. These coordinates are required to be
     * three dimensional.
     */
    static inline const_coordinate_iterator coordsBegin( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*!
     * \brief Return the const iterator to the end of the node coordinate
     * block in this mesh. These coordinates are requried to be three
     * dimensional. 
     */
    static inline const_coordinate_iterator coordsEnd( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }
    //@}


    //@{
    //! Mesh element concepts.
    /*!
     * \brief Return the element type for this mesh (DTK enum).
     */
    static inline std::size_t elementType( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*! 
     * \brief Return the element topology for this mesh (DTK enum).
     */
    static inline std::size_t elementTopology( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*! 
     * \brief Return the number of nodes that constructs an individual element
     * in this mesh.
     */
    static inline std::size_t nodesPerElement( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*! 
     * \brief Return the const iterator to the beginning of the element handle
     * block in this mesh.
     */
    static inline const_handle_iterator elementsBegin( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*! 
     * \brief Return the const iterator to the end of the element handle block
     * in this mesh.
     */
    static inline const_handle_iterator elementsEnd( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*! 
     * \brief Return the const iterator to the beginning of the element
     * connectivity block in this mesh. 
     */
    static inline const_handle_iterator connectivityBegin( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }

    /*! 
     * \brief Return the const iterator to the end of the element connectivity
     * block in this mesh.
     */
    static inline const_handle_iterator connectivityEnd( const T& mesh )
    { return UndefinedMeshTraits<T>::notDefined(); }
    //@}
};

} // end namespace DataTransferKit

#endif // end DTK_MESHTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshTraits.hpp
//---------------------------------------------------------------------------//
