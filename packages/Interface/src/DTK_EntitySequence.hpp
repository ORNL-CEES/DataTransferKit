//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \brief DTK_EntitySequence.hpp
 * \author Stuart R. Slattery
 * \brief Entity sequence interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYSEQUENCE_HPP
#define DTK_ENTITYSEQUENCE_HPP

#include <iterator>
#include <functional>

#include "DTK_Entity.hpp"

#include <boost/iterator/filter_iterator.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntitySequence
  \brief Entity sequence interface.

  This class provides a mechanism to iterate over a sequence of entities with
  a specified predicate operation for selection.
*/
//---------------------------------------------------------------------------//
class EntitySequence
{
  public:

    //! Iterator typedef.
    typedef std::iterator<std::forward_iterator_tag,Entity>   ForwardIterator;
    typedef std::function<bool(Entity)>                       Predicate;
    typedef boost::filter_iterator<Predicate,ForwardIterator> FilteredForwardIterator;

    /*!
     * \brief Constructor.
     */
    EntitySequence();

    /*!
     * \brief Destructor.
     */
    virtual ~EntitySequence();

    /*!
     * \brief Number of elements in the sequence.
     * \return Number of elements in the sequence.
     */
    virtual std::size_t size() const;

    /*!
     * \brief An iterator assigned to the beginning the sequence.
     * \return An iterator assigned to the beginning the sequence.
     */
    virtual ForwardIterator begin() const;

    /*!
     * \brief An iterator assigned to the end of the sequence.
     * \return An iterator assigned to the end of the sequence.
     */
    virtual ForwardIterator end() const;

    /*!
     * \brief Get an iterator assigned to the beginning of the sequence
     * filtered by a predicate.
     * \param predicate The filtering predicate.
     * \Return A filtered iterator assigned to the beginning of the sequence.
     */
    FilteredForwardIterator 
    filteredBegin( const Predicate& predicate = selectAll ) const;

    /*!
     * \brief Get an iterator assigned to the end of the sequence filtered by
     * a predicate.
     * \param predicate The filtering predicate.
     * \Return A filtered iterator assigned to the end of the sequence.
     */
    FilteredForwardIterator 
    filteredEnd( const Predicate& predicate = selectAll) const;

    /*!
     * \brief Select all entities predicate.
     */
    static inline bool selectAll( Entity entity ) { return true; }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ENTITYSEQUENCE_HPP

//---------------------------------------------------------------------------//
// end DTK_EntitySequence.hpp
//---------------------------------------------------------------------------//
