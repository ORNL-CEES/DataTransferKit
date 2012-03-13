//---------------------------------------------------------------------------//
/*!
 * \file DataTransferKit_CommIndexer.hpp
 * \author Stuart Slattery
 * \brief CommIndexer declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COMMINDEXER_HPP
#define DTK_COMMINDEXER_HPP

#include <map>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//===========================================================================//
/*!
 * \class CommIndexer
 * \brief Map the process ids of a local communicator into a global
 * communicator that encompasses it.
 *
 */
//===========================================================================//

template<class Ordinal>
class CommIndexer
{
  public:

    //@{
    //! Useful typedefs.
    typedef Teuchos::Comm<Ordinal>                         Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>             RCP_Communicator;
    typedef std::map<Ordinal,Ordinal>                      IndexMap;
    //@}

  private:

    // Local to global process id map.
    IndexMap d_l2gmap;

  public:

    // Constructor.
    CommIndexer( RCP_Communicator global_comm, 
		 RCP_Communicator local_comm );

    // Destructor.
    ~CommIndexer();

    // Given a process id in the local communicator, return the distributed
    // object's process id in the global communicator.
    const Ordinal l2g( const Ordinal local_id ) const;

    // Return the size of the local to global map.
    const int size() const
    { return d_l2gmap.size(); }
};

} // end namespace DataTransferKit

#include "DataTransferKit_CommIndexer_Def.hpp"

#endif // end DTK_COMMINDEXER_HPP

//---------------------------------------------------------------------------//
// end DataTransferKit_CommIndexer.hpp
//---------------------------------------------------------------------------//
