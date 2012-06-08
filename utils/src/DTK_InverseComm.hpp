//---------------------------------------------------------------------------//
/*!
 * \file DTK_InverseComm.hpp
 * \author Stuart R. Slattery
 * \brief Inverse communication operator declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INVERSECOMM_HPP
#define DTK_INVERSECOMM_HPP

#include <multimap>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

#include <Tpetra_Map.hpp>

namespace DataTransferKit
{

template<typename Ordinal>
class InverseComm
{

  public:

    //@{
    //! Typedefs.
    typedef Ordinal                              ordinal_type;
    typedef Teuchos::Comm<int>                   CommType;
    typedef Teuchos::RCP<const CommType>         RCP_Comm;
    typedef Tpetra::Map<ordinal_type>            TpetraMap;
    typedef Teuchos::RCP<Tpetra_Map>             RCP_TpetraMap;
    //@}

    // Constructor.
    InverseComm( const RCP_Comm& comm )

    // Destructor.
    ~InverseComm();

    // Given a vector of destination procs compute the corresponding Tpetra
    // map for import.
    RCP_TpetraMap createImportMap( const std::multimap<int,Ordinal>& data );

  private:
   
    // Communicator.
    RCP_Comm d_comm;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_InverseComm_def.hpp"

#endif // end DTK_INVERSECOMM_HPP

//---------------------------------------------------------------------------//
// end DTK_InverseComm.hpp
//---------------------------------------------------------------------------//

