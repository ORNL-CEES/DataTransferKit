//---------------------------------------------------------------------------//
/*!
 * \file DTK_Element.hpp
 * \author Stuart R. Slattery
 * \brief General finite element class definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ELEMENT_HPP
#define DTK_ELEMENT_HPP

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

template<typename HandleType=int>
class Element
{
  private:

    // Element topology (enum).
    std::size_type d_topology;

    // Handle.
    HandleType d_handle;

    // Connectivity.
    Teuchos::ArrayRCP<HandleType> d_connectivity;

  public:

    //! Default constructor.
    Element()
    { /* ... */ }

    //! Connectivity constructor.
    template<class HandleField>
    Element( const int topology,
	     const HandleType handle, 
	     const HandleField &connectivity )
	: d_topology( topology )
	, d_handle( handle )
	, d_connectivity( connectivity )
    { /* ... */ }

    //! Copy constructor.
    Element( const Element<HandleType>& elem );

    //! Copy constructor.
    Element<HandleType>&
    operator=( const Element<HandleType>& elem );

    //! Get the topology.
    int getTopology() const
    { return d_topology; }

    //! Get the handle.
    HandleType getHandle() const 
    { return d_handle; }

    //! Set the handle.
    void setHandle( const HandleType& handle )
    { d_handle = handle; }

    //! Get the connectivity.
    const Teuchos::ArrayRCP<HandleType>& getConnectivity() const
    { return d_connectivity; }

    //! Get the number of nodes constructing this element.
    int getNumNodes() const
    { return d_connectivity.size(); }
};

} // end namespace DataTransferKit

#endif // end DTK_ELEMENT_HPP

//---------------------------------------------------------------------------//
// end DTK_Element.hpp
//---------------------------------------------------------------------------//

