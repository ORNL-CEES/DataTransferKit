module DataTransferKit
  use iso_c_binding

  enum, bind(C)
    enumerator :: DTK_BLOCKED=1, DTK_INTERLEAVED=2
  end enum

  interface

    function DTK_Map_create(comm, src_coord, src_num, src_layout, &
        tgt_coord, tgt_num, tgt_layout, space_dim, options) result(dtk_map) &
        bind(C, name="DTK_Map_create_f")
      use iso_c_binding
      implicit none
      type(c_ptr) :: dtk_map
      integer(kind=c_int), value :: comm
      type(c_ptr), value :: src_coord
      integer(kind=c_size_t), value :: src_num
      integer(kind=c_int), value :: src_layout
      type(c_ptr), value :: tgt_coord
      integer(kind=c_size_t), value :: tgt_num
      integer(kind=c_int), value :: tgt_layout
      integer(kind=c_int), value :: space_dim
      character(kind=c_char) :: options(*)
    end function DTK_Map_create

    subroutine DTK_Map_apply(dtk_map, src_field, src_layout, tgt_field, &
        tgt_layout, field_dim, apply_transpose) &
        bind(C, name="DTK_Map_apply")
      use iso_c_binding
      implicit none
      type(c_ptr), value :: dtk_map
      type(c_ptr), value :: src_field
      integer(kind=c_int), value :: src_layout
      type(c_ptr), value :: tgt_field
      integer(kind=c_int), value :: tgt_layout
      integer(kind=c_int), value :: field_dim
      logical(kind=c_bool), value :: apply_transpose
    end subroutine DTK_Map_apply

    subroutine DTK_Map_delete(dtk_map) &
        bind(C, name="DTK_Map_delete")
      use iso_c_binding
      implicit none
      type(c_ptr), value :: dtk_map
    end subroutine DTK_Map_delete

  end interface

end module
