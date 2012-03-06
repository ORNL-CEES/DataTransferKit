//---------------------------------------------------------------------------//
/*!
 * \file Coupler_CommIndexer.hpp
 * \author Stuart Slattery
 * \brief CommIndexer declaration.
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_COMMINDEXER_HPP
#define COUPLER_COMMINDEXER_HPP

#include <map>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace Coupler
{

//===========================================================================//
/*!
 * \class CommIndexer
 * \brief Map a distributed object's process id in a local to communicator to
 * its process id in a global communicator. This is effectively computing the
 * intersection of the local communicator and global communicator. 
 *
 */
//===========================================================================//

template<class T, class Ordinal=int>
class CommIndexer
{
  public:

    //@{
    //! Useful typedefs.
    typedef Teuchos::RCP<T>                                RCP_T;
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
		 RCP_Communicator local_comm,
		 RCP_T t_ptr );

    // Destructor.
    ~CommIndexer();

    // Given a process id in the local communicator, return the distributed
    // object's process id in the global communicator.
    const Ordinal l2g( const Ordinal local_id ) const;

    // Return the size of the local to global map.
    const int size() const
    { return d_l2gmap.size(); }
};

} // end namespace Coupler

#include "Coupler_CommIndexer_Def.hpp"

#endif // end COUPLER_COMMINDEXER_HPP

//---------------------------------------------------------------------------//
// end Coupler_CommIndexer.hpp
//---------------------------------------------------------------------------//
