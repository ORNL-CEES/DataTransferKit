//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Packing_Utils.i.hh
 * \author Thomas M. Evans
 * \date   Thu Jan  3 11:22:29 2008
 * \brief  Member definitions of class Packing_Utils.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: Packing_Utils.i.hh,v 1.1 2008/01/03 18:26:44 9te Exp $
//---------------------------------------------------------------------------//

#ifndef utils_Packing_Utils_i_hh
#define utils_Packing_Utils_i_hh

namespace denovo
{

//---------------------------------------------------------------------------//
// INLINE AND TEMPLATE MEMBERS FOR PACKER
//---------------------------------------------------------------------------//
/*!
 * \brief Set an allocated buffer to write data into.
 *
 * If compute_buffer_size_mode() is on, this function turns it off.
 *
 * This function accepts an allocated char* buffer.  It assigns begin and end
 * pointers and a mutable position pointer that acts like an iterator.  The
 * Packer will write POD (Plain Old Data) data into this buffer starting at
 * the beginning address of the buffer.  This function must be called before
 * any Packer::pack() calls can be made.
 *
 * Once Packer::set_buffer() is called, all subsequent calls to Packer::pack()
 * will write data incrementally into the buffer set by set_buffer().  To write
 * data into a different buffer, call Packer::set_buffer() again; at this point
 * the Packer no longer has any knowledge about the old buffer.
 *
 * Note, the buffer must be allocated large enough to hold all the data that
 * the client intends to load into it.  There is no memory allocation
 * performed by the Packer class; thus, the buffer cannot be increased in
 * size if a value is written past the end of the buffer.  Optionally, the
 * required buffer size may also be computed using the
 * compute_buffer_size_mode().  See the Packer::pack() function for more
 * details.
 *
 * \param size_in size of the buffer
 * \param buffer  pointer to the char * buffer
 *
 */
void Packer::set_buffer(unsigned int size_in, pointer buffer)
{
    Require (buffer);

    size_mode = false;
    
    // set the size, begin and end pointers, and iterator
    stream_size  = size_in;
    ptr          = &buffer[0];
    begin_ptr    = &buffer[0];
    end_ptr      = begin_ptr + stream_size;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Depending on mode, pack data into a buffer, or compute increment
 * to buffer size.
 *
 * This function's behavior depends on whether in compute_buffer_size_mode(),
 * or not.
 *
 * In compute_buffer_size_mode(), the sizeof(T) operator is used to
 * add the size of the data to the total stream size.  Once this function is
 * called for all of the data to be packed, the size() member function may be
 * used to retrieve the buffer size required.
 *
 * Note that using compute_buffer_size_mode() is optional.  See examples
 * below.
 *
 * Regardless, once the user allocates the buffer, set_buffer() may then be
 * called, which turns off compute_buffer_size_mode (if on).  A call to
 * pack() then actually packs its argument into the buffer.  It also advances
 * the pointer (iterator) location to the next location automatically.  It
 * uses the sizeof(T) operator to get the size of the data; thus, only data
 * where sizeof() has meaning will be properly written to the buffer.
 *
 * Packer::pack() does bounds checking to ensure that the buffer and buffer
 * size defined by Packer::set_buffer() are consistent.  This bounds-checking
 * is always on as the Packer is not normally used in compute-intensive
 * calculations.
 *
 * \param value data of type T to pack into the buffer; the data size must be
 * accessible using the sizeof() operator.
 *
 * Example using compute_buffer_size_mode():
   \code
       double d1 = 5.0, d2 = 10.3;         // data to be packed
       Packer p;
       p.compute_buffer_size_mode();
       p << d1 << d2;                      // computes required size
       vector<char> buffer(p.size());      // allocate buffer
       p.set_buffer(p.size(), &buffer[0]);
       p << d1 << d2;                      // packs d1 and d2 into buffer
   \endcode

 * Example not using compute_buffer_size_mode():
   \code
   double d1 = 5.0, d2 = 10.3;
   Packer p;
   unsigned int bsize = 2 * sizeof(double);  // compute buffer size
   vector<char> buffer(bsize);
   p.set_buffer(bsize, &buffer[0]);
   p << d1 << d2;                            // packs d1 and d2 into buffer
   \endcode
 
 */
template<class T>
void Packer::pack(const T &value)
{
    if ( size_mode )
        stream_size += sizeof(T);
    else
    {
        Require (begin_ptr);
        Ensure (ptr >= begin_ptr);
        Ensure (ptr + sizeof(T) <= end_ptr);
	
        // copy value into the buffer
        std::memcpy(ptr, &value, sizeof(T));
	
        // advance the iterator pointer to the next location
        ptr += sizeof(T);
    }
}

//---------------------------------------------------------------------------//
// INLINE AND TEMPLATE MEMBERS FOR UNPACKER
//---------------------------------------------------------------------------//
/*!
 * \brief Set an allocated buffer to read data from.
 *
 * This function accepts an allocated char* buffer.  It assigns begin and end
 * pointers and a mutable position pointer that acts like an iterator.  The
 * Unpacker will read POD data from this buffer starting at the beginning
 * address of the buffer.  This function must be called before any
 * Unpacker::unpack() calls can be made.
 *
 * Once Unpacker::set_buffer() is called, all subsequent calls to
 * Unpacker::unpack() will read data incrementally from the buffer set by
 * set_buffer().  To read data from a different buffer, call
 * Unpacker::set_buffer() again; at this point the Unpacker no longer has any
 * knowledge about the old buffer.
 *
 * Note, there is no memory allocation performed by the Unpacker class.  Also,
 * the client must know how much data to read from the stream (of course
 * checks can be made telling where the end of the stream is located using
 * the Unpacker::get_ptr(), Unpacker::begin(), and Unpacker::end() functions).
 *
 * \param size_in size of the buffer
 * \param buffer  const_pointer to the char * buffer
 *
 */
void Unpacker::set_buffer(unsigned int size_in, const_pointer buffer)
{
    Require (buffer);
    
    // set the size, begin and end pointers, and iterator
    stream_size  = size_in;
    ptr          = &buffer[0];
    begin_ptr    = &buffer[0];
    end_ptr      = begin_ptr + stream_size;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack data from the buffer.

 * This function unpacks a piece of data (single datum) from the buffer set
 * by Unpacker::set_buffer().  It advances the pointer (iterator) location to
 * the next location automatically.  It uses the sizeof(T) operator to get
 * the size of the data; thus, only data where sizeof() has meaning will be
 * properly read from the buffer.
 *
 * Unpacker::unpack() does bounds checking to ensure that the buffer and
 * buffer size defined by Unpacker::set_buffer() are consistent.  This
 * bounds-checking is always on as this should not be used in computation
 * intensive parts of the code.
 *
 * \param value data of type T to unpack from the buffer; the data size must
 *              be accessible using the sizeof() operator
 */
template<class T>
void Unpacker::unpack(T &value)
{
    Require (begin_ptr);
    Ensure  (ptr >= begin_ptr);
    Ensure  (ptr + sizeof(T) <= end_ptr);

    // copy data into the value reference
    std::memcpy(&value, ptr, sizeof(T));

    // advance the iterator pointer to the next location
    ptr += sizeof(T);
}

} // end namespace denovo

#endif // utils_Packing_Utils_i_hh

//---------------------------------------------------------------------------//
//              end of utils/Packing_Utils.i.hh
//---------------------------------------------------------------------------//
