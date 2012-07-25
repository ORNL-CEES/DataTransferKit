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

//---------------------------------------------------------------------------//
/*!
 * \class CommIndexer
 * \brief Map the process ids of a local communicator into a global
 * communicator that encompasses it.
 */
//---------------------------------------------------------------------------//
class CommIndexer
{
  public:

    //@{
    //! Useful typedefs.
    typedef Teuchos::Comm<int>                             CommType;
    typedef Teuchos::RCP<const CommType>                   RCP_Comm;
    typedef std::map<int,int>                              IndexMap;
    //@}

  private:

    // Local to global process id map.
    IndexMap d_l2gmap;

  public:

    // Constructor.
    CommIndexer( RCP_Comm global_comm, RCP_Comm local_comm );

    // Destructor.
    ~CommIndexer();

    // Given a process id in the local communicator, return the distributed
    // object's process id in the global communicator.
    const int l2g( const int local_id ) const;

    //! Return the size of the local to global map.
    const int size() const
    { return d_l2gmap.size(); }
};

} // end namespace DataTransferKit

#endif // end DTK_COMMINDEXER_HPP

//---------------------------------------------------------------------------//
// end DataTransferKit_CommIndexer.hpp
//---------------------------------------------------------------------------//
