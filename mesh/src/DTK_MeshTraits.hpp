//---------------------------------------------------------------------------//
/*!
 * \brief DTK_MeshTraits.hpp
 * \author Stuart R. Slattery
 * \brief Declaration of mesh traits.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTRAITS_HPP
#define DTK_MESHTRAITS_HPP

#include <iterator>

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
class MeshTraits
{
  public:

    //@{
    //! Typedefs.
    //! Typedef for mesh type.
    typedef T mesh_type;

    //! Typedef for handle type. This type must implement
    //! Teuchos::OrdinalTraits.
    typedef typename T::handle_type handle_type;

    //! Typedef for random access const iterator to node handle values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const handle_type>  
    const_node_iterator;

    //! Typedef for random access const iterator to coordinate values. This
    //! is enforcing a coordinate type of double.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const double>  
    const_coordinate_iterator;

    //! Typedef for random access const iterator to element handle values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const handle_type>  
    const_element_iterator;

    //! Typedef for random access const iterator to connectivity values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const handle_type>  
    const_connectivity_iterator;
    //@}


    //@{
    //! Mesh node concepts.
    /*!
     * \brief Return the const iterator to the beginning of the node handle
     * block in this mesh. 
    */
    static inline const_node_iterator nodesBegin( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*!
     * \brief Return the const iterator to the end of the node handle block in
     * this mesh.
    */ 
    static inline const_node_iterator nodesEnd( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*!
     * \brief Return the const iterator to the beginning of the node
     * coordinate block in this mesh. These coordinates are required to be
     * three dimensional and blocked.
     * ( x0, x1, x2, ... , xN, y0, y1, y2, ... , yN, z0, z1, z2, ... , zN )
     */
    static inline const_coordinate_iterator coordsBegin( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*!
     * \brief Return the const iterator to the end of the node coordinate
     * block in this mesh. These coordinates are requried to be three
     * dimensional and blocked.
     * ( x0, x1, x2, ... , xN, y0, y1, y2, ... , yN, z0, z1, z2, ... , zN )
     */
    static inline const_coordinate_iterator coordsEnd( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }
    //@}


    //@{
    //! Mesh element concepts.
    /*!
     * \brief Return the element type for this mesh (DTK enum).
     */
    static inline std::size_t elementType( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*! 
     * \brief Return the element topology for this mesh (DTK enum).
     */
    static inline std::size_t elementTopology( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0;}

    /*! 
     * \brief Return the number of nodes that constructs an individual element
     * in this mesh. All elements in the mesh must be constructed with the
     * same number of nodes.
     */
    static inline std::size_t nodesPerElement( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the beginning of the element handle
     * block in this mesh.
     */
    static inline const_element_iterator elementsBegin( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the end of the element handle block
     * in this mesh.
     */
    static inline const_element_iterator elementsEnd( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the beginning of the element
     * connectivity block in this mesh. The connectivity entries are required
     * to be blocked. 
     * ( element0( c0 ), element1( c0 ), ... , elementN( c0 ), element0( c1 ),
     * element1( c1 ), ... , elementN( c1 ), ... , elementN( cn ) )
     */
    static inline const_connectivity_iterator connectivityBegin( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the end of the element connectivity
     * block in this mesh. The connectivity entries are required to be blocked. 
     * ( element0( c0 ), element1( c0 ), ... , elementN( c0 ), element0( c1 ),
     * element1( c1 ), ... , elementN( c1 ), ... , elementN( cn ) )
     */
    static inline const_connectivity_iterator connectivityEnd( const T& mesh )
    { UndefinedMeshTraits<T>::notDefined(); return 0; }
    //@}
};

} // end namespace DataTransferKit

#endif // end DTK_MESHTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshTraits.hpp
//---------------------------------------------------------------------------//
