//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in history and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of history code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
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
 * \file DTK_BufferDataTraits.hpp
 * \author Stuart R. Slattery
 * \brief Data buffer type traits definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BUFFERDATA_HPP
#define DTK_BUFFERDATA_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class UndefinedBufferDataTraits
 * \brief Class for undefined data buffer traits. 
 *
 * Will throw a compile-time error if these traits are not specialized.
 */
template<class T>
struct UndefinedBufferDataTraits
{
    static inline void notDefined()
    {
	return T::this_type_is_missing_a_specialization();
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class BufferDataTraits
 * \brief Traits for packing data into buffers.
 */
template<class T>
class BufferDataTraits
{
  public:

    //@{
    //! Typedefs.
    typedef T                                      data_type;
    typedef typename T::ordinal_type               ordinal_type;
    //@}

    /*!
     * \brief Create a data from a buffer.
     */
    static Teuchos::RCP<T> 
    createFromBuffer( const Teuchos::ArrayView<char>& buffer )
    { 
	UndefinedBufferDataTraits<T>::notDefined(); 
	return Teuchos::null;
    }

    /*!
     * \brief Pack the data into a buffer.
     */
    static Teuchos::Array<char> pack( const T& data )
    {
	UndefinedBufferDataTraits<T>::notDefined(); 
	return Teuchos::Array<char>(0);
    }

    /*!
     * \brief Set the byte size of the packed data state.
     */
    static void setByteSize()
    {
	UndefinedBufferDataTraits<T>::notDefined();
    }

    /*!
     * \brief Get the number of bytes in the packed data state.
     */
    static std::size_t getPackedBytes()
    {
	UndefinedBufferDataTraits<T>::notDefined();
	return 0;
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_BUFFERDATA_HPP

//---------------------------------------------------------------------------//
// end DTK_BufferDataTraits.hpp
//---------------------------------------------------------------------------//

