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
C BELONGS IN MPC
C ******************************************************************** C
      function form_tmod_vec(n_sh, beam_id) result(tmod)
        ! Returns the tangent modulus for each node.
        ! @input n_sh: Number of shell elements on interface.
        ! @input beam_id: Beam node tag for the interface.
        ! @returns: The tangent modulus averaged for each node.
        integer, intent(in)   ::  n_sh, beam_id
        integer               ::  asize
        parameter(asize=1000)
        integer               ::  TMOD_ARR_ID, ID_ADJ
        parameter(TMOD_ARR_ID=1, ID_ADJ=3)
        real(8)               ::  tmod_arr(asize)
        real(8)               ::  tmod(n_sh), nelem
        integer               ::  info_arr(3*n_sh), interf_id, ii,
     1                            elem_indexs(3)
        pointer(p_tmod, tmod_arr)
        pointer(p_info, info_arr)
#include <SMAAspUserSubroutines.hdr>
        ! Function start
        interf_id = beam_id + ID_ADJ
        p_tmod = SMAFloatArrayAccess(TMOD_ARR_ID)
        p_info = SMAIntArrayAccess(interf_id)
        tmod = 0.d0
        do ii = 1, n_sh
          elem_indexs(1:3) = info_arr(3*ii-2:3*ii)
          nelem = 0.
          do j = 1, 3
            if (elem_indexs(j) /= 0) then
              tmod(ii) = tmod(ii) + tmod_arr(elem_indexs(j))
              nelem = nelem + 1.
            end if
          end do
          tmod(ii) = tmod(ii) / nelem
        end do
      end function
C ******************************************************************** C

C ******************************************************************** C
C BELONGS IN UMAT
C ******************************************************************** C
      function set_tmod(noel) result(res)
        ! Sets the tangent modulus for the shell elem on the interface.
        ! @input cep: Tangent moduli matrix (plane stress).
        ! @input noel: Element number of the integration point.
        !real                ::  cep(:, :)
        integer             ::  asize
        parameter(asize=1000)
        integer             ::  TMOD_ARR_ID, ID_ARR_ID, DIR_ARR_ID
        parameter(TMOD_ARR_ID=1,ID_ARR_ID=2,DIR_ARR_ID=3)
        real(8)             ::  tmod_arr(asize)
        integer             ::  id_arr(asize),direc_arr(asize),arr_loc,
     1                          narr, direc, res
        pointer(p_tmod, tmod_arr)
        pointer(p_id, id_arr)
        pointer(p_direc, direc_arr)
#include <SMAAspUserSubroutines.hdr>
        ! Function start
        narr = SMAIntArraySize(ID_ARR_ID)
        p_id = SMAIntArrayAccess(ID_ARR_ID)
        arr_loc = loc_interf_elem(noel, id_arr, narr)
        if (arr_loc /= -1) then
          p_direc = SMAIntArrayAccess(DIR_ARR_ID)
          direc = direc_arr(arr_loc)
          p_tmod = SMAFloatArrayAccess(TMOD_ARR_ID)
        ! todo: this only works for one integration point, what about sections in the S4R elem?
          tmod_arr(arr_loc) = ddsdde(direc, direc)
        end if
        res = 0
      end function
C ******************************************************************** C        
      function loc_interf_elem(noel, id_arr, narr) result(res)
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
      end function
C ******************************************************************** C

