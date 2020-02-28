
  function solve_weights(mrow, ncol, amat, bvec) result(wvec)
    ! Returns the weight vector for equilibrium
    ! @input mrow: Number of rows in A.
    ! @input ncol: Number of columns in A.
    ! @input amat: (mrow,ncol) The equilibrium matrix A.
    ! @input bvec: (mrow,1) Right-hand side vector.
    ! @returns: Weights that satisfy equilibrium.
    ! Input and output
    implicit none
    integer, intent(in) ::  mrow, ncol
    real(8), intent(in) ::  amat(mrow, ncol), bvec(mrow)
    real(8)             ::  wvec(ncol)
    ! Internal
    integer             ::  lwork, info, lda, ldb, nrhs
    integer             ::  lwmax
    parameter(lwmax=1000)
    real(8)             ::  work(lwmax)
    external dgels
    wvec = 0.d0
    wvec(1:mrow) = bvec
    lda = mrow
    ldb = ncol
    nrhs = 1
    lwork = -1
    call dgels('N', mrow, ncol, nrhs, amat, lda, wvec, ldb, work, lwork, info)
    lwork = min(lwmax, int(work(1)))
    call dgels('N', mrow, ncol, nrhs, amat, lda, wvec, ldb, work, lwork, info)
  end function
  
! ******************************************************************** !

  function rigid_constraint(n_sh, xdef, xtrans, xcent, normal_vec, psi) result(fvec)
    ! Returns the rigid restraint forces
    ! Input and output
    implicit none
    integer, intent(in) ::  n_sh
    real(8), intent(in) ::  xdef(3, n_sh), xtrans(3, n_sh), xcent(3), normal_vec(3), psi(n_sh)
    real(8)             ::  fvec(3*n_sh)
    ! Internal
    integer             ::  neqn
    parameter(neqn=7)
    real(8)             ::  amat(neqn, 3*n_sh), fdiag(3,3), wvec(3*n_sh), one_vec(3*n_sh)
    integer             ::  ii
    ! Function start
    one_vec = 1.d0
    fdiag = 0.d0
    amat = 0.d0
    do ii = 1, n_sh
      fvec(3*ii-2:3*ii) = xdef(:, ii) - xtrans(:, ii)
      amat(1:3, 3*ii-2:3*ii) = 0.d0
      amat(1, 3*ii-2) = fvec(3*ii-2)
      amat(2, 3*ii-1) = fvec(3*ii-1)
      amat(3, 3*ii) = fvec(3*ii)
      fdiag(1, 1) = fvec(3*ii-2)
      fdiag(2, 2) = fvec(3*ii-1)
      fdiag(3, 3) = fvec(3*ii)
      ! subtract the centroid location to take the moment about the centroid
      amat(4:6, 3*ii-2:3*ii) = matmul(-skew_sym(xdef(1:3, ii) - xcent), fdiag)
      amat(7, 3*ii-2:3*ii) = psi(ii) * dot_product(normal_vec, fvec(3*ii-2:3*ii)) * normal_vec
    end do
    wvec = solve_weights(neqn, 3*n_sh, amat, matmul(amat, one_vec))
    wvec = wvec - one_vec
    do ii = 1, 3*n_sh
      fvec(ii) = fvec(ii) * wvec(ii)
    end do
  end function

! ******************************************************************** !


  function weights_def_config(n_sh, u_shell, uc, weights_ref, xdef, xc, normal_vec, psi) result(weights_def)
    ! Returns the weights that are in equilibrium in the deformed configuation.
    ! @param n_sh: Number of shell nodes on interface.
    ! @param u_shell: Displacement of the shell nodes
    ! @param uc: Centroid displacement.
    ! @param weights_ref: Weights in the reference configuration.
    ! @param xdef: Shell nodes in the deformed configuration.
    ! @param xc: Position of the centroid in the deformed configuration.
    ! @param normal_vec: Unit vector normal to the cross-section interface in the deformed config.
    ! @param psi: Warping fuction for each node.
    ! @return: Weights adjusted for equilbirium in the deformed configuration.
    implicit none
    ! Input and output
    integer, intent(in) ::  n_sh
    real(8), intent(in) ::  u_shell(3, n_sh), uc(3), weights_ref(3 * n_sh, 3), xdef(3, n_sh), xc(3), normal_vec(3), psi(n_sh)
    real(8)             ::  weights_def(3 * n_sh, 3)
    ! Internal
    integer             ::  neqn
    parameter(neqn=8)
    ! parameter(neqn=7)
    real(8)             ::  b_target(neqn), one_vec(3 * n_sh), b(neqn), equil_mat(neqn, 3 * n_sh), &
                            weight_adjustment(3 * n_sh), diag_mat(3, 3)
    integer             ::  i, j
    ! Function definition
    diag_mat = 0.d0
    one_vec = 1.d0
    equil_mat = 0.d0
    do i = 1, 3
      do j = 1, n_sh
          ! Force equilibrium conditions
          equil_mat(1, 3*j-2) = weights_ref(3*j-2, i)
          equil_mat(2, 3*j-1) = weights_ref(3*j-1, i)
          equil_mat(3, 3*j)   = weights_ref(3*j,   i)
          ! Moment equilibrium conditions
          diag_mat(1, 1) = weights_ref(3*j-2, i)
          diag_mat(2, 2) = weights_ref(3*j-1, i)
          diag_mat(3, 3) = weights_ref(3*j,   i)
          equil_mat(4:6, 3*j-2:3*j) = matmul(-skew_sym(xdef(1:3, j) - xc), diag_mat)
          ! Bimoment equilibrium
          equil_mat(7, 3*j-2:3*j) = psi(j) * dot_product(normal_vec, weights_ref(3*j-2:3*j, i)) * normal_vec
          ! Centroid displacement consitency
          equil_mat(8, 3*j-2:3*j) = weights_ref(3*j-2:3*j, i) * u_shell(1:3, j)
      end do
      b_target = 0.d0
      b_target(i) = 1.d0
      b_target(8) = uc(i)
      b = b_target - matmul(equil_mat, one_vec)
      
      weight_adjustment = one_vec + solve_weights(neqn, 3 * n_sh, equil_mat, b)
      print *, 'weight adj. norm, i = ', i, norm2(weight_adjustment)
      do j = 1, n_sh
          weights_def(3*j-2, i) = weights_ref(3*j-2, i) * weight_adjustment(3*j-2)
          weights_def(3*j-1, i) = weights_ref(3*j-1, i) * weight_adjustment(3*j-1)
          weights_def(3*j  , i) = weights_ref(3*j,   i) * weight_adjustment(3*j)
      end do
    end do
  end function

! ******************************************************************** !
