! Best-fit beam-to-shell MPC coupling for wide-flange cross-sections
!
! Couples beam and shell element domains for wide-flange cross-sections 
! using an MPC approach. 
!
! Notes:
!     - The beam element must have 7 DOF (3 disp, 3 rot, 1 warping)
!     - The coupling is only defined for wide-flange (I-shaped) cross-
!     sections
!     - The centroid of the shell domain must be initially aligned with
!     the beam node
!     - The beam 1 axis is the "strong axis", the "2" is the weak axis
!     - Intel MKL __MUST__ be linked with Abaqus for this subroutine to
!     work
!
! References:
!     [1] Hartloper, Lignos and de Sousa (2019), Best-fit Beam-to-shell 
!     Coupling for Wide-flange Cross-sections
!     [2] Mostafa and Sivaselvan (2014), On best-fit corotated frames 
!     for 3D continuum finite elements
!     [3] Horn (1987), Closed-form solution of absolute orientation 
!     using unit quaternions
!
! Written by: A Hartloper, EPFL, alexander.hartloper@epfl.ch

! Intel pre-processor command for free-form input in .for files
!DIR$ FREEFORM

! Include (with the pre-processor) the modules/subroutines used with MPC
#include 'mpc_modules.f90'
#include 'umat_subr.f90'
#include 'uexternaldb_subr.f90'
!DIR$ FREEFORM

subroutine mpc(ue, a, jdof, mdof, n, jtype, x, u, uinit, maxdof, lmpc, kstep, kinc, time, &
               nt, nf, temp, field, ltran, tran)
  use mpc_modules
!DIR$ NOFREEFORM
      include 'ABA_PARAM.INC'
!DIR$ FREEFORM
  
  ! Parameter definitions
  dimension ue(mdof), a(mdof, mdof, n), jdof(mdof, n), x(6, n), &
    u(maxdof, n), uinit(maxdof, n), time(2), temp(nt, n), &
    field(nf, nt, n), ltran(n), tran(3, 3, n)
  ! Internal variables
  ! Defines the spatial DOF and size of quaternions
  integer               :: ndim, quatdim
  parameter                (ndim = 3, quatdim = 4)
  ! Number of shell nodes
  integer               :: n_shell
  real(8)               :: n_reciprocal
  ! Loop counters
  integer               :: i, i_sh
  ! Nodal weight factors
  real(8)               :: total_area, a_node
  ! Warping amplitude
  real(8)               :: w_amp
  ! Centroids (deformed and initial)
  real(8)               :: xc(ndim), xc0(ndim), u_cent(ndim)
  ! Rotation matrices (to get deformed, reference)
  real(8)               :: r_mat(ndim, ndim), r_ref(ndim, ndim), r_tot(ndim, ndim)
  ! Orientation normal to beam plane
  real(8)               :: t_def(ndim)
  ! Arrays that depend on the number of nodes
  real(8), allocatable  :: c_mat(:, :), d_mat(:, :), g_mat(:, :), &
                           q_mat(:, :), x_shell_def(:, :), &
                           x_shell(:, :), u_shell(:, :), &
                           w_lin(:), psi_all(:), weights(:, :), &
                           x_shell_ref(:, :), q_mat2(:, :)
  ! For the rotation quaternion and computations
  real(8)               ::  b_mat(quatdim, quatdim), &
                            beta(quatdim, quatdim), &
                            beta_T(quatdim, quatdim)
  real(8)               :: op_quat(quatdim)
  real(8)               :: lambda
  real(8)               :: lam_q(quatdim+1)
  ! Linearized displacment
  integer               ::  k_inner, n_inner
  real(8)               ::  tol_inner, uc_prev(3)
  real(8), allocatable  :: disp_lin(:, :),bend_tang(:, :), &
                            my_torque(:, :), rotw(:, :), &
                            wtmod(:), tmod(:), tmod0(:), xy_bar(:, :), f_nodes(:, :)
  real(8)               ::  temp33(3, 3), wtmod_tot, r1, ksec(4, 4), usec(4, 4)
  ! Defines the max number of inner iterations to do and the tolerance
  parameter                 (n_inner=10, tol_inner=1.d-8)
  ! Section properties
  real(8)               ::  section_props(4), d_cl
  real(8)               ::  rigid_scale
  real(8)               ::  sec_props(4), sp_dimensional(4)
  ! Numerical parameters
  real(8)               :: one, two, zero
  parameter                (one=1.0d0, two=2.0d0, zero=0.d0)
  integer               ::  stall_var, nothing_var, interf_id
  ! For the cross-section properties
  pointer(p_secprops, section_props)
!DIR$ NOFREEFORM
#include <SMAAspUserSubroutines.hdr>
!DIR$ FREEFORM
  
  ! Allocate variable size arrays
  n_shell = n - 1
  allocate(x_shell(ndim, n_shell))
  allocate(u_shell(ndim, n_shell))
  allocate(x_shell_def(ndim, n_shell))
  allocate(x_shell_ref(ndim, n_shell))
  allocate(c_mat(ndim, n_shell))
  allocate(d_mat(ndim, n_shell))
  allocate(g_mat(ndim, ndim * n_shell))
  allocate(q_mat(ndim, ndim * n_shell))
  allocate(psi_all(n_shell))
  allocate(w_lin(ndim * n_shell))
  allocate(weights(5, n_shell))
  allocate(disp_lin(ndim, ndim * n_shell))
  allocate(bend_tang(ndim, n_shell))
  allocate(my_torque(ndim, n_shell))
  allocate(tmod(n_shell))
  allocate(tmod0(n_shell))
  allocate(wtmod(n_shell))
  allocate(rotw(4, n_shell))
  allocate(q_mat2(ndim, ndim * n_shell))
  allocate(xy_bar(2, n_shell))
  allocate(f_nodes(4, n_shell))
  allocate(disp_lin2(ndim, ndim * n_shell))
  allocate(w_lin2(ndim * n_shell))

! ******************************************************************** !
! Subroutine definition
! ******************************************************************** !
  
  
  print *, 'kinc', kinc
  
  ! For debugging
  stall_var = 0
  !       print *, 'time', time(1)
  !       if (abs(time(1)-0.57) < 0.0001) then
  !         stall_var = 1
  !         do while (stall_var == 1)
  !           nothing_var = 0
  !         end do
  !       end if
  
  ! Retrieve the cross-section properties from the global array
  ! section_props = [d, bf, tf, tw]
  p_secprops = SMAFloatArrayAccess(2 * jtype + 3)
  
  ! Calculate the centroid locations
  ! The first node is the beam node, the shell nodes follow
  x_shell = x(1:3, 2:n)
  u_shell = u(1:3, 2:n)
  x_shell_def = x_shell + u_shell
  xc0(:) = zero
  do i = 1, n_shell
    xc0(:) = xc0(:) + x_shell(:, i)
  end do
  n_reciprocal = one / n_shell
  xc0 = xc0 * n_reciprocal

  ! Find the reference orientation
  do i = 1, n_shell
    c_mat(:, i) = x_shell(:, i) - xc0
  end do
  r_ref = ini_config(c_mat, n_shell)
  ! Compute the warping function at each node
  psi_all = warp_fun(c_mat, n_shell, r_ref)

  ! Get the tangent
  !tmod = form_tmod_vec(n_shell, jtype)
  tmod = 1.d0
  tmod0 = 1.d0

  ! Calculate the nodal weights for displacement linearization
  ! todo: pass x_shell_ref where needed
  x_shell_ref = matmul(transpose(r_ref), c_mat)
  weights = calc_weights(x_shell_ref, n_shell, section_props, tmod0)
  total_area = zero
  wtmod_tot = 0.d0
  do i = 1, n_shell
    wtmod(i) = weights(1, i) * tmod(i)
    total_area = total_area + weights(1, i)
    wtmod_tot = wtmod_tot + wtmod(i)
  end do  
  weights(1, :) = wtmod(:)
  total_area = wtmod_tot

  ! Initialize the linearization to the arithmetic mean
  temp33 = zero
  temp33(1, 1) = n_reciprocal
  temp33(2, 2) = n_reciprocal
  temp33(3, 3) = n_reciprocal
  do i = 1, n_shell
    disp_lin(1:3, 3*i-2:3*i) = temp33
  end do
  
  ! Iterate until convergence using the displacement linearization for the centroid
  u_cent(1:3) = zero
  do k_inner = 1, n_inner
    ! Compute the centroid location
    uc_prev = u_cent
    u_cent(1:3) = zero
    do i = 1, n_shell
      ! todo: if the deformed config does not affect the centroid location, do we need to bother with the extra computations here?
      ! if (k_inner == 1) then
      !   ! Want to use mean of points for the first inner iteration
      !   temp33 = disp_lin(1:3, 3*i-2:3*i)
      !   u_cent(1:3) = u_cent(1:3) + matmul(temp33, u_shell(1:3, i))
      ! else
      !   u_cent(1) = u_cent(1) + dot_product(weights_def(3*i-2:3*i, 1), u_shell(1:3, i))
      !   u_cent(2) = u_cent(2) + dot_product(weights_def(3*i-2:3*i, 2), u_shell(1:3, i))
      !   u_cent(3) = u_cent(3) + dot_product(weights_def(3*i-2:3*i, 3), u_shell(1:3, i))
      ! end if
      temp33 = disp_lin(1:3, 3*i-2:3*i)
      u_cent(1:3) = u_cent(1:3) + matmul(temp33, u_shell(1:3, i))
    end do
    
    xc = xc0 + u_cent
    ! Center the deformed config on the centroid
    do i = 1, n_shell
      d_mat(:, i) = x_shell_def(:, i) - xc
    end do

    ! Compute the optimal rotation quaternion
    b_mat(:, :) = zero
    do i = 1, n_shell
      beta=left_quat(zero,d_mat(:,i)) - right_quat(zero,c_mat(:,i))
      beta_T = transpose(beta)
      b_mat = b_mat + wtmod(i) / wtmod_tot * matmul(beta_T,beta)
    end do
    lam_q = calc_opquat(b_mat)
    
    lambda = lam_q(1)
    op_quat = lam_q(2:5)
    ! Make the first entry always positive
    op_quat = sign(one, op_quat(1)) * op_quat  
    r_mat = rot_mat(op_quat)

    ! Compute the linearized centroid displacement
    disp_lin = calc_disp_lin(n_shell, weights, total_area, r_ref, r_mat)
    
    ! Check for convergene of the inner iterations
    if (norm2(u_cent-uc_prev)/(norm2(uc_prev)+1.d-9)<tol_inner) then
      exit
    end if
  end do  ! over k_inner
    
  print *, 'u'
  print *, u_shell
  print *, 'ref 2'
  print *, disp_lin(2, :)
  print *, 'r tot'
  print *, transpose(r_tot)
  print *, 'u cent'
  print *, u_cent
  print *, 'x def'
  print *, x_shell_def
  print *, 'def 2'
  print *, weights_def(:, 2)
  
  ! Compute the linearized rotation matrix
  g_mat = calc_g(op_quat, n_shell, b_mat, c_mat, lambda, r_mat)
  q_mat = calc_q(op_quat, n_shell, g_mat)

  ! Compute the warping amplitude
  t_def = matmul(r_mat, r_ref(:, 3))
  w_amp = warp_amp(x_shell,x_shell_def,u_cent,t_def,r_mat,psi_all)

  ! Compute the linearized warping 
  w_lin = calc_lin_w(x_shell_ref,n_shell,r_ref(:,3),psi_all,r_mat, weights(1, :))

  ! Assign the A submatrix for the beam node and set active DOF
  forall(i = 1:maxdof) a(i, i, 1) = one
  forall(i = 1:maxdof) jdof(i, 1) = i
  
  ! todo: clean this up with actual values (if actually used)
  ! bend_tang = tang_force_b(n_shell, x_shell_ref, r_ref, r_mat, 25.0d0, 25.0d0, 28.d0, 14.5d0, 300.d0, 472.d0)

  ! Compute the inelastic nodal forces based on the normal stresses
  xy_bar = normalized_coords(n_shell, x_shell_ref(1:2, :), section_props(1), section_props(2))
  sec_props = compute_section_props(n_shell, xy_bar, weights(1, :))
  sp_dimensional = compute_section_props(n_shell, x_shell_ref(1:2, :), weights(1, :))
  ! ksec = section_tangent(n_shell, xy_bar, weights(1, :), tmod, sec_props)
  ksec = section_tangent(n_shell, xy_bar, weights(1, :), tmod0, sec_props)
  usec = calc_section_deform(ksec, sp_dimensional, section_props(1), section_props(2))
  ! f_nodes = nodal_forces(n_shell, xy_bar, weights(1, :), tmod, usec)
  f_nodes = nodal_forces(n_shell, xy_bar, weights(1, :), tmod0, usec)

  ! Rotation linearization (inelastic)
  ! Torque weights
  d_cl = section_props(1) - section_props(3)
  my_torque = torquer(n_shell, x_shell_ref, r_mat, weights(1, :), d_cl)
  ! todo: use the torque weights
  ! at the moment replace the Mx and My with the ones computed using the tangent modulus
  rotw = rot_weights(n_shell, x_shell_ref, weights(1, :), tmod0)
  rotw(1:2, :) = f_nodes(2:3, :)
  q_mat2 = lin_rot(n_shell, r_ref, r_mat, rotw)

  ! Rotate the inelastic nodal forces for displacement
  ! todo: Fix this hack that replaces the area with the tangent normalized one
  weights(1, :) = f_nodes(1, :)
  disp_lin2 = calc_disp_lin(n_shell, weights, 1.d0, r_ref, r_mat)

  ! Linearized warping
  do i = 1, n_shell
    w_lin2(3*i-2:3*i) = t_def * f_nodes(4, i)
  end do
  
  ! Assign the A submatrices for the shell nodes
  do i = 2, n
    ! i_sh corresponds to the shell nodes since i starts at 2
    i_sh = i - 1
    ! a_node = weights(1, i_sh) / total_area
    a_node = wtmod(i_sh) / wtmod_tot
    ! Displacement constraints
    !a(1:3, 1:3, i) = -disp_lin(1:3, 3*i_sh-2:3*i_sh)
    a(1:3, 1:3, i) = -weights_def(3*i_sh-2:3*i_sh, 1:3)
    ! Rotation constraints
    q_mat(1:3, 3*i_sh-2:3*i_sh)=q_mat(1:3, 3*i_sh-2:3*i_sh)*a_node
    a(4:6, 1:3, i) = -q_mat(1:3, 3*i_sh-2:3*i_sh)
    !a(4:6, 1:3, i) = -q_mat2(1:3, 3*i_sh-2:3*i_sh)

    ! Address tangential due to bending moment
    !a(4, 1:3, i) = a(4, 1:3, i) + bend_tang(1:3, i_sh)
    ! Replace the torque
    !a(6, 1:3, i) = my_torque(1:3, i_sh)

    ! Warping constraint
    a(7, 1:3, i) = -w_lin(3*i_sh-2:3*i_sh)
    ! Set the active DOF (displacement) in the shell elements
    jdof(1, i) = 1
    jdof(2, i) = 2
    jdof(3, i) = 3
  end do
  
  ! Update the beam node DOF exactly
  ue(1:3) = u_cent
  ue(4:6) = extract_rotation(op_quat)
  ue(7) = w_amp
  
  ! if (x(3,1) > 2500.) then
  ! print *, 'ucent'
  ! print *, u_cent
  ! print *, 'rot'
  ! print *, extract_rotation(op_quat)
  ! print *, 'disp_lin'
  ! print *, disp_lin(2, :)
  ! end if
  
  return

end subroutine
