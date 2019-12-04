C ******************************************************************** C
C ******************************************************************** C
      subroutine uexternaldb(lop, lrestart, time, dtime, kstep, kinc)
        ! Sets-up the global arrays for the tangent modulus.
        include 'ABA_PARAM.INC'
        dimension time(2)
        ! Start subroutine
        integer ::  read_unit, read_result
        parameter(read_unit=99)
#include <SMAAspUserSubroutines.hdr>
        if (lop == 0) then
          open(unit=read_unit, file='D:/Abaqus_Files/MPC_Coupling
     1_HEB500_4m/run_files/interface_props.txt', status='old')
          read_result = input_file_reader(read_unit)
          close(unit=read_unit)
        end if
        return
      contains
C ******************************************************************** C
      function input_file_reader(read_unit) result(res)
        ! Reads input file and sets-up the global arrays.
        ! @input read_unit: Unit to read from.
        !
        ! Notes:
        !   - Only allows for up to 1000 entries in any global array
        !   - This max value can be modified by changing asize
        integer, intent(in)   ::  read_unit
        integer               ::  res
        integer               ::  n_val, ios, global_id, j
        character(5)          ::  keyword, elem_key, interf_key, 
     1                            end_file_key
        integer               ::  TMOD_ARR_ID, ID_ARR_ID, DIR_ARR_ID
        parameter(TMOD_ARR_ID=1,ID_ARR_ID=2,DIR_ARR_ID=3)
        integer, allocatable  ::  elem_array(:, :), interf_array(:, :)
        integer               ::  arr_size
        parameter(asize=1000)
        integer               ::  id_arr(asize), dir_arr(asize), 
     1                            interf_arr(asize)
        real(8)               ::  tmod_arr(asize)
        pointer(ptr_tmod, tmod_arr)
        pointer(ptr_id, id_arr)
        pointer(ptr_dir, dir_arr)
        pointer(ptr_iarr, interf_arr)
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
            ptr_tmod = SMAFloatArrayCreate(TMOD_ARR_ID, n_val, 1.0)
            ptr_id = SMAIntArrayCreate(ID_ARR_ID, n_val, -1)
            ptr_dir = SMAIntArrayCreate(DIR_ARR_ID, n_val, -1)
            do j = 1, n_val
              id_arr(j) = elem_array(1, j)
              dir_arr(j) = elem_array(2, j)
            end do
          
          else if (keyword == interf_key) then
            allocate(interf_array(3, n_val))
            read(read_unit, *) global_id
            interf_array = interface_reader(read_unit, n_val)
            ! Create the global array
            interf_id = global_id + 3
            ptr_iarr = SMAIntArrayCreate(interf_id, 3 * n_val, -1)
            do j = 1, n_val
              interf_arr(3*j-2:3*j) = interf_array(1:3, j)
            end do
            deallocate(interf_array)
            
          else if (keyword == end_file_key) then
            exit
          end if
          
          if (ios /= 0) exit
        end do
        res = 0
      end function
C ******************************************************************** C
      function elements_reader(iunit, n_elem) result(elem_id_dir)
        ! Reads and returns an array of the element id and direction.
        ! @param iunit: Open unit to read from.
        ! @param n_elem: Number of entries to read.
        ! @returns: (2,n_elem) Array with id and direction.
        integer, intent(in) ::  iunit, n_elem
        integer             ::  elem_id_dir(2, n_elem)
        integer             ::  i
        ! Function start
        do i = 1, n_elem
          read(iunit, *) elem_id_dir(1:2, i)
        end do
      end function
C ******************************************************************** C
      function interface_reader(iunit, n_node) result(n2e_map)
        ! Reads and returns an array of the node-to-element map.
        ! @param iunit: Open unit to read from.
        ! @param n_node: Number of entries to read.
        ! @returns: (3,n_node) Array with 3 indices for each node.
        integer, intent(in) ::  iunit, n_node
        integer             ::  n2e_map(3, n_node)
        integer             ::  i
        ! Function start
        do i = 1, n_node
          read(iunit, *) n2e_map(1:3, i)
        end do
      end function
      end   ! SUBROUTINE UEXTERNALDB
C ******************************************************************** C
C ******************************************************************** C


      pure function form_tmod_vec(n_sh, info_arr_id) result(tmod)
        ! Returns the tangent modulus for each node.
        ! @input n_sh: Number of shell elements on interface.
        ! @input info_arr_id: ID of the global array assoc. with
        !   coupling.
        ! @returns: The tangent modulus averaged for each node.
        integer, intent(in)   ::  info_arr_id
        real(8)               ::  tmod_arr(:)
        real(8)               ::  tmod(n_sh), nelem
        integer               ::  info_arr(3, n_sh), TMOD_ARR_ID, n_tm
        parameter(TMOD_ARR_ID=1)
        pointer(p_tmod, tmod_arr)
        pointer(p_info, info_arr)
        ! Function start
        ! n_tm = SMAIntArraySize(TMOD_ARR_ID)
        ! allocate(tmod_arr(n_tm))
        p_tmod = SMAFloatArrayAccess(TMOD_ARR_ID)
        p_info = SMAIntArrayAccess(info_arr_id)

        tmod = 0.d0
        do i = 1, n_sh
          elem_indexs = info_arr(:, i)
          nelem = 0.
          do j = 1, 3
            if (elem_indexs(j) .not. 0) then
              tmod(i) = tmod(i) + tmod_arr(elem_indexs(j))
              nelem = nelem + 1.
            end if
          end do
          tmod(i) = tmod(i) / nelem
        end do
      end

      pure function set_tmod(arr_loc, dir, cep)
        ! Sets the tangent modulus for the shell elem on the interface.
        ! @input arr_loc: Index for the global tangent modulus array.
        ! @input dir: Direction to use from the tangent moduli matrix.
        ! @input cep: Tangent moduli matrix (plane stress).
        integer, intent(in) ::  info_arr_id
        real(8), intent(in) ::  cep(3, 3)
        real(8)             ::  tmod_arr(:)
        integer             ::  dir_arr(:), TMOD_ARR_ID
        pointer(p_tmod, tmod_arr)
        ! todo: put these in some import?
        parameter(TMOD_ARR_ID=1)
        ! Function start
        p_tmod = SMAFloatArrayAccess(TMOD_ARR_ID)
        tmod_arr(arr_loc) = cep(dir, dir)
      end


        
      pure function loc_interf_elem(noel, id_arr, narr) result(res)
        ! Returns the index in id_arr that matches noel.
        ! @input noel: ID of element to find in id_arr.
        ! @input id_arr: Array to search.
        ! @input narr: Length of id_arr.
        ! @returns: i such that id_arr(i)==noel, or -1 if not found.
        integer, intent(in) ::  narr, noel, id_arr(narr)
        integer             ::  res, ii
        res = -1
        do ii = 1, narr
          ! todo: improve the linear search if we have a sorted array
          if (id_arr(ii) == noel) then
            res = ii
            exit
          end if
        end do
      end

