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
! TODOs:
!   - Clean-up the variables that are defined / used
!   - Look at efficiency, joining loops together in functions
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
  real(8), allocatable  :: c_mat(:, :), d_mat(:, :), &
                           rot_lin(:, :), x_shell_def(:, :), &
                           x_shell(:, :), u_shell(:, :), &
                           warp_lin(:), psi_all(:), weights(:, :), &
                           x_shell_ref(:, :), node_areas(:), shear_weights(:, :)
  integer, allocatable  ::  node_classes(:)
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
  real(8)               ::  section_props(4), d_cl, bf
  real(8)               ::  rigid_scale, mesh_sizes(2)
  real(8)               ::  sec_props(4), area_and_moi(4)
  ! Numerical parameters
  real(8)               :: one, two, zero
  parameter                (one=1.0d0, two=2.0d0, zero=0.d0)
  integer               ::  stall_var, nothing_var, interf_id
  ! For the cross-section properties global arrays
  pointer(p_secprops, section_props)
!DIR$ NOFREEFORM
#include <SMAAspUserSubroutines.hdr>
!DIR$ FREEFORM
  
  ! Allocate variable size arrays
  ! todo: remove all the entries that aren't needed
  n_shell = n - 1
  allocate(x_shell(ndim, n_shell))
  allocate(u_shell(ndim, n_shell))
  allocate(x_shell_def(ndim, n_shell))
  allocate(x_shell_ref(ndim, n_shell))
  allocate(c_mat(ndim, n_shell))
  allocate(d_mat(ndim, n_shell))
  allocate(rot_lin(ndim, ndim * n_shell))
  allocate(psi_all(n_shell))
  allocate(warp_lin(ndim * n_shell))
  allocate(weights(4, n_shell))
  allocate(disp_lin(ndim, ndim * n_shell))
  allocate(tmod(n_shell))
  allocate(rotw(4, n_shell))
  allocate(xy_bar(2, n_shell))
  allocate(f_nodes(4, n_shell))
  allocate(node_areas(n_shell))
  allocate(shear_weights(n_shell, 3))
  allocate(node_classes(n_shell))

! ******************************************************************** !
! Subroutine definition
! ******************************************************************** !
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
  x_shell_ref = matmul(transpose(r_ref), c_mat)

  ! Compute section properties
  d_cl = section_props(1) - section_props(3)
  bf = section_props(2)
  node_classes = classify_nodes(n_shell, x_shell_ref, d_cl, bf)
  mesh_sizes = compute_mesh_size(n_shell, x_shell_ref, section_props, node_classes)
  node_areas = compute_node_areas(n_shell, x_shell_ref, section_props, node_classes, mesh_sizes)
  xy_bar = normalized_coords(n_shell, x_shell_ref(1:2, :), section_props(1), section_props(2))
  area_and_moi = compute_section_props(n_shell, x_shell_ref(1:2, :), node_areas)
  total_area = area_and_moi(1)
  psi_all = warp_fun(c_mat, n_shell, r_ref)
  ! The tangent modulus is assumed to be elastic at all shell nodes
  tmod = 1.d0
  
  ! Calculate the elastic shear weights
  shear_weights = shear_factors(x_shell_ref, n_shell, section_props, mesh_sizes, node_classes)
  ! Assume elastic, ksec is identity, else can compute the section stiffness
  forall(i = 1:4) ksec(i, i) = one
  usec = calc_section_deform(ksec, area_and_moi, section_props(1), section_props(2))
  f_nodes = nodal_forces(n_shell, xy_bar, node_areas, tmod, usec)
  weights(1, :) = node_areas
  weights(2:4, :) = transpose(shear_weights)

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
      b_mat = b_mat + node_areas(i) / total_area * matmul(beta_T,beta)
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
 
  ! Rotation linearization using nodal force method
  rotw(1:2, :) = f_nodes(2:3, :)
  rotw(3:4, :) = torquer(n_shell, x_shell_ref, node_areas, d_cl)
  rot_lin = lin_rot(n_shell, r_ref, r_mat, rotw)
 
  ! Compute the warping amplitude
  t_def = matmul(r_mat, r_ref(:, 3))
  w_amp = warp_amp(x_shell, x_shell_def, u_cent, t_def, r_mat, psi_all)
  ! Linearized warping
  do i = 1, n_shell
    warp_lin(3*i-2:3*i) = t_def * f_nodes(4, i)
  end do
  
  ! Assign the A submatrix for the beam node and set active DOF
  forall(i = 1:maxdof) a(i, i, 1) = one
  forall(i = 1:maxdof) jdof(i, 1) = i
  
  ! Assign the A submatrices for the shell nodes
  do i = 2, n
    ! i_sh corresponds to the shell nodes since i starts at 2
    i_sh = i - 1
    ! Displacement constraints
    a(1:3, 1:3, i) = -disp_lin(1:3, 3*i_sh-2:3*i_sh)
    ! Rotation constraints
    a(4:6, 1:3, i) = -rot_lin(1:3, 3*i_sh-2:3*i_sh)
    ! Warping constraint
    a(7, 1:3, i) = -warp_lin(3*i_sh-2:3*i_sh)
    ! Set the active DOF (displacement) in the shell elements
    jdof(1, i) = 1
    jdof(2, i) = 2
    jdof(3, i) = 3
  end do
  
  ! Update the beam node DOF exactly
  ue(1:3) = u_cent
  ue(4:6) = extract_rotation(op_quat)
  ue(7) = w_amp
  return

end subroutine
