module utc
  use, intrinsic :: ISO_C_BINDING
  implicit none

  ! PUBLIC METHODS AND TYPES
  public :: UserTestClass
  public :: check_registry

  ! TYPES
  type, bind(C) :: UserTestClass
      integer(c_int), public :: space_dim
      integer(c_size_t), public :: size_1
      integer(c_size_t), public :: size_2
      integer(c_int), public :: offset
      type(c_ptr), public :: field_name
      type(c_ptr), public :: data
  end type

  private
  interface
   function swigc_check_registry(farg1, farg2) &
       bind(C, name="swigc_check_registry") &
       result(fresult)
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT) :: fresult
     character(C_CHAR) :: farg1
     type(C_PTR), value :: farg2
   end function
  end interface

contains
  function check_registry(name, handle) &
     result(fresult)
   use, intrinsic :: ISO_C_BINDING
   integer(C_INT) :: fresult
   character(kind=C_CHAR, len=*) :: name
   type(C_PTR) :: handle
   fresult = swigc_check_registry(name, handle)
  end function

end module utc

module c_f_string_module
  implicit none
  private
  public :: c_f_string
  public :: f_c_string_ptr
contains
  function c_f_string(c_str) result(f_str)
    use, intrinsic :: iso_c_binding
    type(c_ptr), intent(in) :: c_str
    character(kind=c_char, len=:), pointer :: f_str
    character(kind=c_char), dimension(:), pointer :: chars
    integer(c_size_t) :: length, i
    interface
      ! steal std c library function rather than writing our own.
      function strlen(s) bind(c, name='strlen')
        use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
        implicit none
        type(c_ptr), intent(in), value :: s
        integer(c_size_t) :: strlen
      end function strlen
    end interface
    length = strlen(c_str)
    call c_f_pointer(c_str, chars, [length])
    allocate(character(kind=c_char, len=length) :: f_str)
    do i = 1, length
      f_str(i:i) = chars(i)
    enddo
  end function c_f_string
  subroutine f_c_string_ptr(f_str, c_str)
    use, intrinsic :: iso_c_binding
    character(kind=c_char, len=256) :: f_str
    type(c_ptr), value :: c_str

    integer :: i
    character(kind=c_char), dimension(:), pointer :: chars

    call c_f_pointer(c_str, chars, [256])
    do i = 1, len(f_str)
      chars(i) = f_str(i:i)
    enddo
  end subroutine f_c_string_ptr

end module c_f_string_module

module x
  use, intrinsic :: ISO_C_BINDING
  use datatransferkit
  use utc
  use c_f_string_module

  implicit none

contains
  !---------------------------------------------------------------------------
  ! User functions.
  !---------------------------------------------------------------------------
  ! Get the size parameters for building a node list.
  subroutine node_list_size(user_data, space_dim, local_num_nodes) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_int), intent(out) :: space_dim
    integer(c_size_t), intent(out) :: local_num_nodes

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    space_dim = u%space_dim
    local_num_nodes = u%size_1
  end subroutine

  !---------------------------------------------------------------------------
  ! Get the data for a node list.
  subroutine node_list_data(user_data, coordinates) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    real(c_double), dimension(*), intent(out) :: coordinates

    integer :: i, j
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      do j = 1, u%space_dim
        coordinates((j-1) * u%size_1 + i) = (i-1) + (j-1) + u%offset
      end do
    end do
  end subroutine

  !---------------------------------------------------------------------------
  ! Get the size parameters for building a bounding volume list.
  subroutine bounding_volume_list_size(user_data, space_dim, local_num_volumes) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_int), intent(out) :: space_dim
    integer(c_size_t), intent(out) :: local_num_volumes

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    space_dim = u%space_dim
    local_num_volumes = u%size_1
  end subroutine

  !---------------------------------------------------------------------------
  ! Get the data for a bounding volume list.
  subroutine bounding_volume_list_data(user_data, bounding_volumes) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    real(c_double), dimension(*), intent(out) :: bounding_volumes

    integer :: i, j, h, index
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      do j = 1, u%space_dim
        do h = 1, 2
          index = u%size_1 * u%space_dim * (h-1) + u%size_1 * (j-1) + i
          bounding_volumes(index) = (i-1) + (j-1) + (h-1) + u%offset
        end do
      end do
    end do
  end subroutine

  !---------------------------------------------------------------------------
  ! Get the size parameters for building a polyhedron list.
  subroutine polyhedron_list_size(user_data, space_dim, local_num_nodes, &
      local_num_faces, total_nodes_per_face, local_num_cells, total_faces_per_cell) &
      BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_int), intent(out) :: space_dim
    integer(c_size_t), intent(out) :: local_num_nodes
    integer(c_size_t), intent(out) :: local_num_faces
    integer(c_size_t), intent(out) :: total_nodes_per_face
    integer(c_size_t), intent(out) :: local_num_cells
    integer(c_size_t), intent(out) :: total_faces_per_cell

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    space_dim = u%space_dim
    local_num_nodes = u%size_1
    local_num_faces = u%size_1
    total_nodes_per_face = u%size_1
    local_num_cells = u%size_1
    total_faces_per_cell = u%size_1
  end subroutine

  !---------------------------------------------------------------------------
  ! Get the data for a polyhedron list.
  ! FIXME: there is no c_unsigned
  subroutine polyhedron_list_data(user_data, coordinates, faces, nodes_per_face, &
      cells, faces_per_cell, face_orientation) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    real(c_double), dimension(*), intent(out) :: coordinates
    integer(c_int), dimension(*), intent(out) :: faces
    integer(c_int), dimension(*), intent(out) :: nodes_per_face
    integer(c_int), dimension(*), intent(out) :: cells
    integer(c_int), dimension(*), intent(out) :: faces_per_cell
    integer(c_int), dimension(*), intent(out) :: face_orientation

    integer :: i, j
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      do j = 1, u%space_dim
        coordinates((j-1) * u%size_1 + i) = (i-1) + (j-1) + u%offset
      end do
      faces(i) = (i-1) + u%offset
      nodes_per_face(i) = (i-1) + u%offset
      cells(i) = (i-1) + u%offset
      faces_per_cell(i) = (i-1) + u%offset
      face_orientation(i) = 1
    end do
  end subroutine

  !---------------------------------------------------------------------------
  ! Get the size parameters for building a cell list.
  ! FIXME: there is no c_unsigned
  subroutine cell_list_size(user_data, space_dim, local_num_nodes, &
      local_num_cells, total_cell_nodes) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_int), intent(out) :: space_dim
    integer(c_size_t), intent(out) :: local_num_nodes
    integer(c_size_t), intent(out) :: local_num_cells
    integer(c_size_t), intent(out) :: total_cell_nodes

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    space_dim = u%space_dim
    local_num_nodes = u%size_1
    local_num_cells = u%size_1
    total_cell_nodes = u%size_1
  end subroutine

  !---------------------------------------------------------------------------
  ! Get the data for a cell list.
  subroutine cell_list_data(user_data, coordinates, cells, cell_topologies) &
      BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    real(c_double), dimension(*), intent(out) :: coordinates
    integer(c_int), dimension(*), intent(out) :: cells
    integer(kind(DTK_CellTopology)), dimension(*), intent(out) :: cell_topologies

    integer :: i, j
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      do j = 1, u%space_dim
        coordinates((j-1) * u%size_1 + i) = (i-1) + (j-1) + u%offset
      end do
      cells(i) = (i-1) + u%offset
      cell_topologies(i) = DTK_TET_4
    end do
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the size parameters for a boundary.
  subroutine boundary_size(user_data, local_num_faces) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_size_t), intent(out) :: local_num_faces

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    local_num_faces = u%size_1
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the data for a boundary.
  ! FIXME: there is no c_unsigned
  subroutine boundary_data(user_data, boundary_cells, cell_faces_on_boundary) &
      BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_int), dimension(*), intent(out) :: boundary_cells
    integer(c_int), dimension(*), intent(out) :: cell_faces_on_boundary

    integer :: i
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      cell_faces_on_boundary(i) = (i-1) + u%offset
      boundary_cells(i) = (i-1) + u%offset
    end do
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the size parameters for building a cell list.
  subroutine adjacency_list_size(user_data, total_adjacencies) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_size_t), intent(out) :: total_adjacencies

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    total_adjacencies = u%size_1
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the data for a boundary.
  ! FIXME: there is no c_unsigned
  subroutine adjacency_list_data(user_data, global_cell_ids, &
          adjacent_global_cell_ids, adjacencies_per_cell) &
      BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_long_long), dimension(*), intent(out) :: global_cell_ids
    integer(c_long_long), dimension(*), intent(out) :: adjacent_global_cell_ids
    integer(c_int), dimension(*), intent(out) :: adjacencies_per_cell

    integer :: i
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      global_cell_ids(i) = (i-1) + u%offset
      adjacent_global_cell_ids(i) = i-1
      adjacencies_per_cell(i) = 1
    end do
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the size parameters for a degree-of-freedom id map with a single
  ! number of dofs per object.
  ! FIXME: there is no c_unsigned
  subroutine dof_map_size(user_data, local_num_dofs, local_num_objects, &
      dofs_per_object) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_size_t), intent(out) :: local_num_dofs
    integer(c_size_t), intent(out) :: local_num_objects
    integer(c_int), intent(out) :: dofs_per_object

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    local_num_dofs = u%size_1
    local_num_objects = u%size_1
    dofs_per_object = u%size_2
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the data for a degree-of-freedom id map with a single number of
  ! dofs per object.
  subroutine dof_map_data(user_data, global_dof_ids, object_dof_ids, &
      discretization_type) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_long_long), dimension(*), intent(out) :: global_dof_ids
    integer(c_int), dimension(*), intent(out) :: object_dof_ids
    type(c_ptr), value :: discretization_type
    character(kind=c_char, len=256) :: fstring = "unit_test_discretization"//C_NULL_CHAR

    integer :: i, j
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      global_dof_ids(i) = (i-1) + u%offset
      do j = 1, u%size_2
        object_dof_ids((j-1) * u%size_1 + i) = (i-1) + (j-1) + u%offset
      end do
    end do

    call f_c_string_ptr(fstring, discretization_type)
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the size parameters for a degree-of-freedom id map with a
  ! multiple number of dofs per object (e.g. mixed topology cell lists or
  ! polyhedron lists).
  subroutine mixed_topology_dof_map_size(user_data, local_num_dofs, &
      local_num_objects, total_dofs_per_object) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_size_t), intent(out) :: local_num_dofs
    integer(c_size_t), intent(out) :: local_num_objects
    integer(c_size_t), intent(out) :: total_dofs_per_object

    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    local_num_dofs = u%size_1
    local_num_objects = u%size_1
    total_dofs_per_object = u%size_1
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the data for a multiple object degree-of-freedom id map
  ! (e.g. mixed topology cell lists or polyhedron lists).
  ! FIXME: there is no c_unsigned
  subroutine mixed_topology_dof_map_data(user_data, global_dof_ids, &
      object_dof_ids, dofs_per_object, discretization_type) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    integer(c_long_long), dimension(*), intent(out) :: global_dof_ids
    integer(c_int), dimension(*), intent(out) :: object_dof_ids
    integer(c_int), dimension(*), intent(out) :: dofs_per_object
    type(c_ptr), value :: discretization_type
    character(kind=c_char, len=256) :: fstring = "unit_test_discretization"//C_NULL_CHAR

    integer :: i
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    do i = 1, u%size_1
      global_dof_ids(i) = (i-1) + u%offset
      object_dof_ids(i) = (i-1) + u%offset
      dofs_per_object(i) = u%size_2
    end do
    call f_c_string_ptr(fstring, discretization_type)
  end subroutine
  !---------------------------------------------------------------------------
  ! Get the size parameters for a field. Field must be of size
  ! local_num_dofs in the associated dof_id_map.
  ! FIXME: there is no c_unsigned
  subroutine field_size(user_data, field_name, &
      field_dimension, local_num_dofs) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    type(c_ptr), value, intent(in) :: field_name
    integer(c_int), intent(out) :: field_dimension
    integer(c_size_t), intent(out) :: local_num_dofs

    character(kind=c_char,len=:), pointer :: f_field_name
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    f_field_name => c_f_string(field_name)
    ! Here one could do actions depening on the name, but in the tests we
    ! simply ignore it

    field_dimension = u%space_dim
    local_num_dofs = u%size_1
  end subroutine
  !---------------------------------------------------------------------------//
  ! Pull data from application into a field.
  subroutine pull_field_data(user_data, field_name, &
      field_dofs) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    type(c_ptr), value, intent(in) :: field_name
    real(c_double), dimension(*), intent(out) :: field_dofs

    integer :: i, j
    character(kind=c_char,len=:), pointer :: f_field_name
    type(UserTestClass), pointer :: u
    real(c_double), dimension(:), pointer :: udata

    call c_f_pointer(user_data, u)
    call c_f_pointer(u%data, udata, [u%size_1*u%size_2])

    f_field_name => c_f_string(field_name)
    ! Here one could do actions depening on the name, but in the tests we
    ! simply ignore it

    do i = 1, u%size_1
      do j = 1, u%space_dim
        field_dofs((j-1) * u%size_1 + i) = udata((j-1) * u%size_1 + i)
      end do
    end do
  end subroutine
  !---------------------------------------------------------------------------//
  ! Push data from application into a field.
  subroutine push_field_data(user_data, field_name, &
      field_dofs) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    type(c_ptr), value, intent(in) :: field_name
    real(c_double), dimension(*), intent(in) :: field_dofs

    integer :: i, j
    character(kind=c_char,len=:), pointer :: f_field_name
    type(UserTestClass), pointer :: u
    real(c_double), dimension(:), pointer :: udata

    call c_f_pointer(user_data, u)
    call c_f_pointer(u%data, udata, [u%size_1*u%size_2])

    f_field_name => c_f_string(field_name)
    ! Here one could do actions depening on the name, but in the tests we
    ! simply ignore it

    do i = 1, u%size_1
      do j = 1, u%space_dim
        udata((j-1) * u%size_1 + i) = field_dofs((j-1) * u%size_1 + i)
      end do
    end do
  end subroutine
  !---------------------------------------------------------------------------//
  ! Evaluate a field at a given set of points in a given set of objects.
  subroutine evaluate_field(user_data, field_name, &
      num_points, &
      evaluation_points, object_ids, values) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(c_ptr), value :: user_data
    type(c_ptr), value, intent(in) :: field_name
    integer(c_size_t), value, intent(in) :: num_points
    real(c_double), dimension(*), intent(in) :: evaluation_points
    integer(c_int), dimension(*), intent(in) :: object_ids
    real(c_double), dimension(*), intent(out) :: values

    integer :: i, j
    ! character(kind=c_char,len=:), pointer :: f_field_name
    type(UserTestClass), pointer :: u

    call c_f_pointer(user_data, u)

    ! f_field_name => c_f_string(field_name)
    ! Here one could do actions depening on the name, but in the tests we
    ! simply ignore it

    do i = 1, num_points
      do j = 1, u%space_dim
        values((j-1) * num_points + i) = &
          evaluation_points((j-1) * num_points + i) + object_ids(i)
      end do
    end do
  end subroutine

  !---------------------------------------------------------------------------
  ! Tests
  !---------------------------------------------------------------------------
  function test_node_list(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_NODE_LIST_SIZE_FUNCTION, C_FUNLOC(node_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_NODE_LIST_DATA_FUNCTION, C_FUNLOC(node_list_data), C_LOC(u))

    rv = check_registry( "test_node_list"//C_NULL_CHAR, dtk_handle )
  end function

  function test_bounding_volume_list(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_BOUNDING_VOLUME_LIST_SIZE_FUNCTION, C_FUNLOC(bounding_volume_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_BOUNDING_VOLUME_LIST_DATA_FUNCTION, C_FUNLOC(bounding_volume_list_data), C_LOC(u))

    rv = check_registry( "test_bounding_volume_list"//C_NULL_CHAR, dtk_handle )
  end function

  function test_polyhedron_list(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_POLYHEDRON_LIST_SIZE_FUNCTION, C_FUNLOC(polyhedron_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_POLYHEDRON_LIST_DATA_FUNCTION, C_FUNLOC(polyhedron_list_data), C_LOC(u))

    rv = check_registry( "test_polyhedron_list"//C_NULL_CHAR, dtk_handle )
  end function

  function test_multiple_topology_cell(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_CELL_LIST_SIZE_FUNCTION, C_FUNLOC(cell_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_CELL_LIST_DATA_FUNCTION, C_FUNLOC(cell_list_data), C_LOC(u))

    rv = check_registry( "test_multiple_topology_cell"//C_NULL_CHAR, dtk_handle )
  end function

  function test_boundary(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_BOUNDARY_SIZE_FUNCTION, C_FUNLOC(boundary_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_BOUNDARY_DATA_FUNCTION, C_FUNLOC(boundary_data), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_CELL_LIST_SIZE_FUNCTION, C_FUNLOC(cell_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_CELL_LIST_DATA_FUNCTION, C_FUNLOC(cell_list_data), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_POLYHEDRON_LIST_SIZE_FUNCTION, C_FUNLOC(polyhedron_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_POLYHEDRON_LIST_DATA_FUNCTION, C_FUNLOC(polyhedron_list_data), C_LOC(u))

    rv = check_registry( "test_boundary"//C_NULL_CHAR, dtk_handle )
  end function

  function test_adjacency_list(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_ADJACENCY_LIST_SIZE_FUNCTION, C_FUNLOC(adjacency_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_ADJACENCY_LIST_DATA_FUNCTION, C_FUNLOC(adjacency_list_data), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_CELL_LIST_SIZE_FUNCTION, C_FUNLOC(cell_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_CELL_LIST_DATA_FUNCTION, C_FUNLOC(cell_list_data), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_POLYHEDRON_LIST_SIZE_FUNCTION, C_FUNLOC(polyhedron_list_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_POLYHEDRON_LIST_DATA_FUNCTION, C_FUNLOC(polyhedron_list_data), C_LOC(u))

    rv = check_registry( "test_adjacency_list"//C_NULL_CHAR, dtk_handle )
  end function

  function test_single_topology_dof(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_DOF_MAP_SIZE_FUNCTION, C_FUNLOC(dof_map_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_DOF_MAP_DATA_FUNCTION, C_FUNLOC(dof_map_data), C_LOC(u))

    rv = check_registry( "test_single_topology_dof"//C_NULL_CHAR, dtk_handle )
  end function

  function test_multiple_topology_dof(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION, &
      C_FUNLOC(mixed_topology_dof_map_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION, &
      C_FUNLOC(mixed_topology_dof_map_data), C_LOC(u))

    rv = check_registry( "test_multiple_topology_dof"//C_NULL_CHAR, dtk_handle )
  end function

  function test_field_push_pull(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_FIELD_SIZE_FUNCTION, C_FUNLOC(field_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_PULL_FIELD_DATA_FUNCTION, C_FUNLOC(pull_field_data), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_PUSH_FIELD_DATA_FUNCTION, C_FUNLOC(push_field_data), C_LOC(u))

    rv = check_registry( "test_field_push_pull"//C_NULL_CHAR, dtk_handle )
  end function

  function test_field_eval(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_FIELD_SIZE_FUNCTION, C_FUNLOC(field_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_EVALUATE_FIELD_FUNCTION, C_FUNLOC(evaluate_field), C_LOC(u))

    rv = check_registry( "test_field_eval"//C_NULL_CHAR, dtk_handle )
  end function

  function test_missing_function(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_NODE_LIST_SIZE_FUNCTION, C_FUNLOC(node_list_size), C_LOC(u))

    rv = check_registry( "test_missing_function"//C_NULL_CHAR, dtk_handle )
  end function

  function test_too_many_functions(dtk_handle, u) result(rv)
    type(c_ptr), value :: dtk_handle
    type(UserTestClass), target :: u
    integer :: rv

    call DTK_set_user_function( dtk_handle, DTK_DOF_MAP_SIZE_FUNCTION, C_FUNLOC(dof_map_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_DOF_MAP_DATA_FUNCTION, C_FUNLOC(dof_map_data), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_SIZE_FUNCTION, &
      C_FUNLOC(mixed_topology_dof_map_size), C_LOC(u))
    call DTK_set_user_function( dtk_handle, DTK_MIXED_TOPOLOGY_DOF_MAP_DATA_FUNCTION, &
      C_FUNLOC(mixed_topology_dof_map_data), C_LOC(u))

    rv = check_registry( "test_too_many_functions"//C_NULL_CHAR, dtk_handle )
  end function

end module x

program main

#include "DTK_APIConstants.h"

  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING
  use datatransferkit
  use x
  use utc
  use mpi
  implicit none

  integer :: ierr

  integer(c_int) :: my_rank, num_procs
  integer(kind(DTK_MemorySpace)) :: memory_space
  character(256), pointer :: f_field_name
  real(c_double), dimension(:), pointer :: udata

  ! Use target so that later we can call C_LOC on it
  type(UserTestClass), target :: u
  integer :: rv

  integer :: i
  character(len=32) :: opt, optarg

  ! Use multiple handles due to lack of scope in Fortran
  type(c_ptr) :: dtk_handle_1, dtk_handle_2, dtk_handle_3
  type(c_ptr) :: dtk_handle_4, dtk_handle_5, dtk_handle_6
  type(c_ptr) :: dtk_handle_7, dtk_handle_8, dtk_handle_9
  type(c_ptr) :: dtk_handle_10, dtk_handle_11, dtk_handle_12
  type(c_ptr) :: dtk_handle_13, dtk_handle_14

  rv = 0

  u%space_dim = SPACE_DIM
  u%size_1 = SIZE_1
  u%size_2 = SIZE_2
  u%offset = OFFSET

  allocate(f_field_name)
  f_field_name = FIELD_NAME//C_NULL_CHAR
  u%field_name = c_loc(f_field_name)

  allocate(udata(u%size_1 * u%space_dim))
  u%data = c_loc(udata)

  ! Initialize MPI subsystem
  call MPI_INIT(ierr)
  if (ierr /= 0) then
    write(*,*) "MPI failed to init"
    stop 1
  endif

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

  memory_space = DTK_HOST_SPACE

  do i = 1, command_argument_count()
    call get_command_argument(i, opt)

    select case (opt)
    case ('-s', '--space')
      call get_command_argument(i+1, optarg)
      if ( optarg == 'serial' ) then
        memory_space = DTK_HOST_SPACE
      elseif ( optarg == 'openmp' ) then
        memory_space = DTK_HOST_SPACE
      elseif ( optarg == 'cuda' ) then
        memory_space = DTK_CUDAUVM_SPACE
      else
        write(*,*) "Unknown memory space\n"
        stop 1
      endif
    case ('-h', '--help')
      if ( my_rank .eq. 0 ) then
        write(*,*) "Usage: cmdline [-s <serial|openmp|cuda>] [-h]\n"
      endif
      stop
    end select
  end do

  if ( my_rank .eq. 0 ) then
    write(*,*) "DTK version: ", DTK_version()
    write(*,*) "DTK hash: ", DTK_git_commit_hash()
  endif

  call DTK_initialize()
  rv = ior(rv, merge(0, 1, DTK_is_initialized() ))

  !{
    dtk_handle_1 = DTK_create_user_application( memory_space )
    rv = ior(rv, merge(0, 1, DTK_is_valid_user_application( dtk_handle_1 )))
    call DTK_destroy_user_application( dtk_handle_1 )
    rv = ior(rv, merge(1, 0, DTK_is_valid_user_application( dtk_handle_1 )))
  !}
  !{
    rv = ior(rv, merge(1, 0, DTK_is_valid_user_application (dtk_handle_2 )))
  !}
  !{
    dtk_handle_3 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_node_list( dtk_handle_3, u ))
    call DTK_destroy_user_application( dtk_handle_3 )
  !}
  !{
    dtk_handle_4 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_bounding_volume_list( dtk_handle_4, u ))
    call DTK_destroy_user_application( dtk_handle_4 )
  !}
  !{
    dtk_handle_5 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_polyhedron_list( dtk_handle_5, u ))
    call DTK_destroy_user_application( dtk_handle_5 )
  !}
  !{
    dtk_handle_6 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_multiple_topology_cell( dtk_handle_6, u ))
    call DTK_destroy_user_application( dtk_handle_6 )
  !}
  !{
    dtk_handle_7 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_boundary( dtk_handle_7, u ))
    call DTK_destroy_user_application( dtk_handle_7 )
  !}
  !{
    dtk_handle_8 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_adjacency_list( dtk_handle_8, u ))
    call DTK_destroy_user_application( dtk_handle_8 )
  !}
  !{
    dtk_handle_9 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_single_topology_dof( dtk_handle_9, u ))
    call DTK_destroy_user_application( dtk_handle_9 )
  !}
  !{
    dtk_handle_10 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_multiple_topology_dof( dtk_handle_10, u ))
    call DTK_destroy_user_application( dtk_handle_10 )
  !}
  !{
    dtk_handle_11 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_field_push_pull( dtk_handle_11, u ))
    call DTK_destroy_user_application( dtk_handle_11 )
  !}
  !{
    dtk_handle_12 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_field_eval( dtk_handle_12, u ))
    call DTK_destroy_user_application( dtk_handle_12 )
  !}
  !{
    dtk_handle_13 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_missing_function( dtk_handle_13, u ))
    call DTK_destroy_user_application( dtk_handle_13 )
  !}
  !{
    dtk_handle_14 = DTK_create_user_application( memory_space )
    rv = ior(rv, test_too_many_functions( dtk_handle_14, u ))
    call DTK_destroy_user_application( dtk_handle_14 )
  !}

  call DTK_finalize()

  if (my_rank .eq. 0) then
    if (rv .eq. 0) then
      write(*,*) "End Result: TEST PASSED"
    else
      write(*,*) "End Result: TEST FAILED"
    endif
  endif


  deallocate(f_field_name)
  deallocate(udata)

  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)

end program
