//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
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
#ifndef DTK_NONAME_H
#define DTK_NONAME_H

#include "stdbool.h"
#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

// Convenience typedef for an opaque pointer
typedef void DTK_Map;

// Possible data layouts:
// * Blocked
//     X1 X2 ... Xn Y1 Y2 ... Yn
// * Interleaved
//     X1 Y1 X2 Y2 ... Xn Yn
typedef enum data_layout { DTK_BLOCKED=1, DTK_INTERLEAVED=2 } DTK_Data_layout;

//----------------------------------------------------------------------------//
DTK_Map* DTK_Map_create( MPI_Comm        comm,
                         double const*   src_coord,
                         unsigned        src_num,
                         DTK_Data_layout src_layout,
                         double const*   tgt_coord,
                         unsigned        tgt_num,
                         DTK_Data_layout tgt_layout,
                         int             space_dim,
                         char const*     options );

//----------------------------------------------------------------------------//
void DTK_Map_apply( DTK_Map*        dtk_map,
                    double const*   src_field,
                    DTK_Data_layout src_layout,
                    double*         tgt_field,
                    DTK_Data_layout tgt_layout,
                    int             field_dim,
                    bool            transpose );

//----------------------------------------------------------------------------//
void DTK_Map_delete( DTK_Map * dtk_map );

//----------------------------------------------------------------------------//
DTK_Map* DTK_Map_create_f( MPI_Fint        fint,
                           double const*   src_coord,
                           unsigned        src_num,
                           DTK_Data_layout src_layout,
                           double const*   tgt_coord,
                           unsigned        tgt_num,
                           DTK_Data_layout tgt_layout,
                           int             space_dim,
                           char const*     options );


#ifdef __cplusplus
}
#endif

#endif // DTK_NONAME_H
