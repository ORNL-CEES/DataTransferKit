//---------------------------------------------------------------------------//
/*!
 * \file DTK_Map3dConsistentInterpolation.hpp
 * \author Stuart R. Slattery
 * \brief Map declaration for a 3D consistent interpolation scheme.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MAP3DCONSISTENTINTERPOLATION_HPP
#define DTK_MAP3DCONSISTENTINTERPOLATION_HPP

#include "DTK_Map.hpp"
#include <DTK_Rendezvous.hpp>
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

template<class DataSource, class DataTarget>
class Map3dConsistentInterpolation : public Map
{

  public:

    //@{
    //! Typedefs.
    typedef DataSource                    DS;
    typedef DS::node_field_type           source_node_field_type;
    typedef DS::node_type                 source_node_type;
    typedef DS::element_field_type        source_element_field_type;
    typedef DS::element_type              source_element_type;
    typedef DS::data_field_type           source_data_field_type;
    typedef DS::data_type                 source_data_type;
    typedef DataTarget                    DT;
    typedef DT::node_field_type           target_node_field_type;
    typedef DT::node_type                 target_node_type;
    typedef DT::data_field_type           target_data_field_type;
    typedef DT::data_type                 target_data_type;
    typedef Teuchos::Comm<int>            CommType;
    typedef Teuchos::RCP<const CommType>  RCP_Comm;
    //@}

    // Constructor.
    Map3dConsistentInterpolation();

    // Destructor.
    ~Map3dConsistentInterpolation();

    // Setup the map.
    void setup();

    // Apply the map.
    void apply();

  private:


};

} // end namespace DataTransferKit

#endif // end DTK_MAP3DCONSISTENTINTERPOLATION_HPP

//---------------------------------------------------------------------------//
// end DTK_Map3dConsistentInterpolation.hpp
//---------------------------------------------------------------------------//

