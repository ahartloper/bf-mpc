!DIR$ FREEFORM

! Notes:
!   - This subroutine creates global arrays that are used to:
!     1) Store the cross-section dimensions associated with each interface
!     2) Allow for transfer of the tangent modulus computed in UMAT subroutines to the MPC subroutine
!   Functions must exist in the UMAT and MPC subroutines for this communication
!   Use "git show umat_edits_incl" for the tag that indicates the last commit with functions that existed 
!   in the UVCplanestress UMAT subroutine.
!   The arrays that allow for transfer of the information are still created in this subroutine,
!   but not used in the MPC's current form.


subroutine uexternaldb(lop, lrestart, time, dtime, kstep, kinc)
  ! Sets-up the global arrays for the tangent modulus.
  include 'ABA_PARAM.INC'
  dimension time(2)
  ! Start subroutine
  integer ::  read_unit, read_result
  parameter(read_unit=99)
!DIR$ NOFREEFORM
#include <SMAAspUserSubroutines.hdr>
!DIR$ FREEFORM
  if (lop == 0) then
    open(unit=read_unit, file='D:/Abaqus_Files/MPC_Coupling_HEB500_4m/run_files/interface_props.txt', status='old')
    read_result = input_file_reader(read_unit)
    close(unit=read_unit)
  end if
return


contains

! ******************************************************************** !

  function elements_reader(iunit, n_elem) result(elem_id_dir)
    ! Reads and returns an array of the element id and direction.
    ! @param iunit: Open unit to read from.
    ! @param n_elem: Number of entries to read.
    ! @returns: (2,n_elem) Array with id and direction.
    implicit none
    integer, intent(in) ::  iunit, n_elem
    integer             ::  elem_id_dir(2, n_elem)
    integer             ::  i
    ! Function start
    do i = 1, n_elem
      read(iunit, *) elem_id_dir(1:2, i)
    end do
  end function
  
! ******************************************************************** !

  function interface_reader(iunit, n_node) result(n2e_map)
    ! Reads and returns an array of the node-to-element map.
    ! @param iunit: Open unit to read from.
    ! @param n_node: Number of entries to read.
    ! @returns: (3,n_node) Array with 3 indices for each node.
    implicit none
    integer, intent(in) ::  iunit, n_node
    integer             ::  n2e_map(3, n_node)
    integer             ::  i
    ! Function start
    do i = 1, n_node
      read(iunit, *) n2e_map(1:3, i)
    end do
  end function

! ******************************************************************** !

  function input_file_reader(read_unit) result(res)
    ! Reads input file and sets-up the global arrays.
    ! @input read_unit: Unit to read from.
    !
    ! Notes:
    !   - Only allows for up to 1000 entries in any global array
    !   - This max value can be modified by changing asize
    !   - Global arrays:
    !     1: Tangent modulus for each element
    !     2: Element number
    !     3: Element direction
    !     <2*node+2>: Node to element map for interface
    !     <2*node+3>: Cross-section properties for interface
    implicit none
    integer, intent(in)   ::  read_unit
    integer               ::  res
    integer               ::  n_val, ios, global_id, j, interf_node_elem_arr_id, &
                              section_prop_arr_id
    character(5)          ::  keyword, elem_key, interf_key, end_file_key
    integer               ::  TMOD_ARR_ID, ID_ARR_ID, DIR_ARR_ID
    parameter(TMOD_ARR_ID=1,ID_ARR_ID=2,DIR_ARR_ID=3)
    integer, allocatable  ::  elem_array(:, :), interf_array(:, :)
    real(8)               ::  sec_props(4), sec_props2(4)
    integer               ::  asize
    parameter(asize=1000)
    integer               ::  id_arr(asize), dir_arr(asize), interf_arr(asize)
    real(8)               ::  tmod_arr(asize)
    pointer(ptr_tmod, tmod_arr)
    pointer(ptr_id, id_arr)
    pointer(ptr_dir, dir_arr)
    pointer(ptr_iarr, interf_arr)
    pointer(ptr_sparr, sec_props2)
    ! Function start
    elem_key = '*elem'
    interf_key = '*intr'
    end_file_key = '*endf'
    ! Process the file
    do
      read(read_unit, *, iostat=ios) keyword, n_val
      
      if (keyword == elem_key) then
        allocate(elem_array(2, n_val))
        elem_array = elements_reader(read_unit, n_val)
        ! Create the global arrays
        ptr_tmod = SMAFloatArrayCreate(TMOD_ARR_ID, 5 * n_val, 1.0)
        ptr_id = SMAIntArrayCreate(ID_ARR_ID, n_val, -1)
        ptr_dir = SMAIntArrayCreate(DIR_ARR_ID, n_val, -1)
        do j = 1, n_val
          id_arr(j) = elem_array(1, j)
          dir_arr(j) = elem_array(2, j)
        end do
      
      else if (keyword == interf_key) then
        allocate(interf_array(3, n_val))
        read(read_unit, *) global_id, sec_props(1:4)
        interf_array = interface_reader(read_unit, n_val)
        ! Create the global array for the node to element map
        interf_node_elem_arr_id = 2 * global_id + 2
        ptr_iarr = SMAIntArrayCreate(interf_node_elem_arr_id, 3 * n_val, -1)
        print *, 'created array id', interf_node_elem_arr_id
        do j = 1, n_val
          interf_arr(3*j-2:3*j) = interf_array(1:3, j)
        end do
        deallocate(interf_array)
        ! Create the global array for the cross-section properties
        section_prop_arr_id = 2 * global_id + 3
        ptr_sparr = SMAFloatArrayCreate(section_prop_arr_id, 4, 0.)
        print *, 'created sec props id', section_prop_arr_id
        do j = 1, 4
          sec_props2(j) = sec_props(j)
        end do
        
      else if (keyword == end_file_key) then
        exit
      end if
      
      if (ios /= 0) exit
    end do
    res = 0
  end function

! ******************************************************************** !

end subroutine  ! SUBROUTINE UEXTERNALDB
