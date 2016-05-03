#ifndef DTK_NONAME_H
#define DTK_NONAME_H

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
typedef enum data_layout { DTK_BLOCKED, DTK_INTERLEAVED } DTK_Data_layout;

//----------------------------------------------------------------------------//
DTK_Map* DTK_Map_create( MPI_Comm        comm,
                         double const*   src_coord,
                         unsigned        src_num,
                         DTK_Data_layout src_layout,
                         double const*   tgt_coord,
                         unsigned        tgt_num,
                         DTK_Data_layout tgt_layout,
                         int             space_dim,
                         char const*     options = "" );

//----------------------------------------------------------------------------//
void DTK_Map_apply( DTK_Map*        dtk_map,
                    double const*   src_field,
                    DTK_Data_layout src_layout,
                    double*         tgt_field,
                    DTK_Data_layout tgt_layout,
                    int             field_dim );

//----------------------------------------------------------------------------//
void DTK_Map_delete( DTK_Map * dtk_map );

#ifdef __cplusplus
}
#endif

#endif // DTK_NONAME_H
