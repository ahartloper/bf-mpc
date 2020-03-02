C Best-fit beam-to-shell MPC coupling for wide-flange cross-sections
C
C Couples beam and shell element domains for wide-flange cross-sections 
C using an MPC approach. 
C
C Notes:
C     - The beam element must have 7 DOF (3 disp, 3 rot, 1 warping)
C     - The coupling is only defined for wide-flange (I-shaped) cross-
C     sections
C     - The centroid of the shell domain must be initially aligned with
C     the beam node
C     - The beam 1 axis is the "strong axis", the "2" is the weak axis
C     - Intel MKL __MUST__ be linked with Abaqus for this subroutine to
C     work
C
C References:
C     [1] Hartloper, Lignos and de Sousa (2019), Best-fit Beam-to-shell 
C     Coupling for Wide-flange Cross-sections
C     [2] Mostafa and Sivaselvan (2014), On best-fit corotated frames 
C     for 3D continuum finite elements
C     [3] Horn (1987), Closed-form solution of absolute orientation 
C     using unit quaternions
C
C Written by: A Hartloper, EPFL, alexander.hartloper@epfl.ch
C 
      subroutine mpc(ue, a, jdof, mdof, n, jtype, x, u, uinit, maxdof,
     1 lmpc, kstep, kinc, time, nt, nf, temp, field, ltran, tran)
      !
      !include 'aba_param.inc'
      !
      dimension ue(mdof), a(mdof, mdof, n), jdof(mdof, n), x(6, n),
     1 u(maxdof, n), uinit(maxdof, n), time(2), temp(nt, n),
     2 field(nf, nt, n), ltran(n), tran(3, 3, n)
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
      real(8)               :: r_mat(ndim, ndim), r_ref(ndim, ndim)
      ! Orientation normal to beam plane
      real(8)               :: t_def(ndim)
      ! Arrays that depend on the number of nodes
      real(8), allocatable  :: c_mat(:, :), d_mat(:, :), g_mat(:, :), 
     1                         q_mat(:, :), x_shell_def(:, :),
     2                         x_shell(:, :), u_shell(:, :), 
     3                         w_lin(:), psi_all(:), weights(:, :),
     4                         x_shell_ref(:, :), q_mat2(:, :)
      ! For the rotation quaternion and computations
      real(8)               ::  b_mat(quatdim, quatdim), 
     1                          beta(quatdim, quatdim), 
     2                          beta_T(quatdim, quatdim)
      real(8)               :: op_quat(quatdim)
      real(8)               :: lambda
      real(8)               :: lam_q(quatdim+1)
      ! Linearized displacment
      integer               ::  k_inner, n_inner
      real(8)               ::  tol_inner, uc_prev(3)
      real(8), allocatable  ::  disp_lin(:, :), bend_tang(:, :),
     1                          my_torque(:, :), rotw(:, :), tmod(:)
      real(8)               ::  temp33(3, 3), tmod_total
      ! Defines the max number of inner iterations to do and the tolerance
      parameter                 (n_inner=10, tol_inner=1.d-8)
      ! Numerical parameters
      real(8)               :: one, two, zero
      parameter                (one=1.0d0, two=2.0d0, zero=0.d0)
      integer               ::  interf_id, inelastic_interface
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
      allocate(weights(4, n_shell))
      allocate(disp_lin(ndim, ndim * n_shell))
      allocate(bend_tang(ndim, n_shell))
      allocate(my_torque(ndim, n_shell))
      allocate(tmod(n_shell))
      allocate(rotw(4, n_shell))
      allocate(q_mat2(ndim, ndim * n_shell))
C ******************************************************************** C
C Subroutine definition
C ******************************************************************** C
      ! Parameter to select if using inelastic interface
      ! if inelastic_interface == 1, use computed tangent modulus
      ! if inelastic_interface /= 1, use constant tangent modulus
      inelastic_interface = 1
      
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
      
      ! Get the tangent modulus for each node
      ! todo: need to find someway of getting the beam ID in here
      if (x(3,1) < 2000) then
        interf_id = 2245
      else
        interf_id = 2247
      end if
      if (inelastic_interface == 1) then
        tmod = form_tmod_vec(n_shell, interf_id)
      else
        tmod = 1.d0
      end if
      
      ! Calculate the nodal weights for displacement linearization
      ! todo: pass x_shell_ref where needed
      x_shell_ref = matmul(transpose(r_ref), c_mat)
      ! todo: need to get all the tf, tw, h, bf, in the main
      weights = calc_weights(x_shell_ref, n_shell, jtype, tmod)
      total_area = zero
      do i = 1, n_shell
        total_area = total_area + weights(1, i)
      end do      
      
      ! Initialize the linearization to the arithmetic mean
      temp33 = zero
      temp33(1, 1) = n_reciprocal
      temp33(2, 2) = n_reciprocal
      temp33(3, 3) = n_reciprocal
      do i = 1, n_shell
        disp_lin(1:3, 3*i-2:3*i) = temp33
      end do
      ! Iterate until convergence using the displacement linearization 
      ! for the centroid
      u_cent(1:3) = zero
      do k_inner = 1, n_inner
        uc_prev = u_cent
        u_cent(1:3) = zero
        do i = 1, n_shell
          temp33 = disp_lin(1:3, 3*i-2:3*i)
          u_cent(1:3) = u_cent(1:3) + matmul(temp33, u_shell(1:3, i))
        end do
        xc = xc0 + u_cent
        do i = 1, n_shell
          d_mat(:, i) = x_shell_def(:, i) - xc
        end do

        ! Compute the optimal rotation quaternion
        b_mat(:, :) = zero
        do i = 1, n_shell
          beta=left_quat(zero,d_mat(:,i))-right_quat(zero,c_mat(:,i))
          beta_T = transpose(beta)
          b_mat = b_mat+weights(1,i)/total_area*matmul(beta_T,beta)
        end do
        lam_q = calc_opquat(b_mat)
        lambda = lam_q(1)
        op_quat = lam_q(2:5)
        ! Make the first entry always positive
        op_quat = sign(one, op_quat(1)) * op_quat  
        r_mat = rot_mat(op_quat)
        
        ! Compute the linearized centroid displacement
        disp_lin=calc_disp_lin(n_shell,weights,total_area,r_ref,r_mat)
        
        ! Check for convergene of the inner iterations
        if (norm2(u_cent-uc_prev)/(norm2(uc_prev)+1.d-9)<tol_inner) then
          exit
        end if
      end do  ! over k_inner

      ! Compute the linearized rotation matrix
      g_mat = calc_g(op_quat, n_shell, b_mat, c_mat, lambda, r_mat)
      q_mat = calc_q(op_quat, n_shell, g_mat)
      
      ! Compute the warping amplitude
      psi_all = warp_fun(c_mat, n_shell, r_ref)
      t_def = matmul(r_mat, r_ref(:, 3))
      w_amp = warp_amp(x_shell,x_shell_def,u_cent,t_def,r_mat,psi_all)
     
      ! Compute the linearized warping 
      w_lin =calc_lin_w(x_shell_ref,n_shell,r_ref(:,3),psi_all,r_mat,
     1                   weights(1, :))
      
      ! Compute the linearized centroid displacement
      disp_lin=calc_disp_lin(n_shell,weights,total_area,r_ref,r_mat)
     
      ! Assign the A submatrix for the beam node and set active DOF
      forall(i = 1:maxdof) a(i, i, 1) = one
      forall(i = 1:maxdof) jdof(i, 1) = i
      
      bend_tang = tang_force_b(n_shell, x_shell_ref, r_ref, r_mat, 
     1            25.0d0, 25.0d0, 28.d0, 14.5d0, 300.d0, 472.d0)
      my_torque = torquer(n_shell, x_shell_ref, r_mat, weights(1, :), 
     1            472.d0)
      
      rotw = rot_weights(n_shell, x_shell_ref, weights(1, :), tmod)
      q_mat2 = lin_rot(n_shell, r_ref, r_mat, rotw)
     
      ! Assign the A submatrices for the shell nodes
      do i = 2, n
        ! i_sh corresponds to the shell nodes since i starts at 2
        i_sh = i - 1
        a_node = weights(1, i_sh) / total_area
        ! Displacement constraints
        a(1:3, 1:3, i) = -disp_lin(1:3, 3*i_sh-2:3*i_sh)
        ! Rotation constraints
        q_mat(1:3, 3*i_sh-2:3*i_sh)=q_mat(1:3, 3*i_sh-2:3*i_sh)*a_node
        !a(4:6, 1:3, i) = -q_mat(1:3, 3*i_sh-2:3*i_sh)
        a(4:6, 1:3, i) = -q_mat2(1:3, 3*i_sh-2:3*i_sh)
        
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
      return
      
      contains
C ******************************************************************** C
C     Skew symmetric matrix from vector
C ******************************************************************** C
      ! Returns the 3x3 skew symmetric matrix from axial vector v
      ! @ input v: Axial vector of length 3
      ! @ returns vss: 3x3 skew symmetric matrix
      pure function skew_sym(v) result(vss)
      real(8), intent(in) :: v(3)
      real(8)             :: vss(3, 3)
      ! Specify the columns of the skew_symmetric matrix
      vss(:, 1) = [ 0.d0, v(3), -v(2) ]
      vss(:, 2) = [ -v(3), 0.d0, v(1) ]
      vss(:, 3) = [ v(2), -v(1), 0.d0 ]
      end function
C ******************************************************************** C
C     Left quaternion representation from 3 element vector
C ******************************************************************** C
      ! Returns the 4x4 left-side matrix representing {v0, v}
      ! @ input v0: Real part of rotation quaternion
      ! @ input v: Imaginary part of rotation quaternion, length 3
      ! @ return ql: Matrix representing left-hand quaternion 
      !              multiplication
      pure function left_quat(v0, v) result(ql)
      real(8), intent(in) :: v0
      real(8), intent(in) :: v(3)
      real(8)             :: ql(4, 4)
      integer             :: ii
      ! 
      ql(:, 1) = [ v0, v(:) ]
      ql(1, 2:4) = -v(:)
      ql(2:4, 2:4) = skew_sym(v)
      forall(ii = 1:3) ql(1 + ii, 1 + ii) = ql(1 + ii, 1 + ii) + v0
      end function
C ******************************************************************** C
C     Right quaternion representation from 3 element vector
C ******************************************************************** C
      ! Returns the 4x4 right-side matrix representing {v0, v}
      ! @ input v0: Real part of rotation quaternion
      ! @ input v: Imaginary part of rotation quaternion, length 3
      ! @ return qr: Matrix representing right-hand quaternion 
      !              multiplication
      pure function right_quat(v0, v) result(qr)
      real(8), intent(in) :: v0
      real(8), intent(in) :: v(3)
      real(8)             :: qr(4, 4)
      integer             :: ii
      !
      qr(:, 1) = [ v0, v(:) ]
      qr(1, 2:4) = -v(:)
      qr(2:4, 2:4) = -skew_sym(v)   
      forall(ii = 1:3) qr(1 + ii, 1 + ii) = qr(1 + ii, 1 + ii) + v0   
      end function
C ******************************************************************** C
C     Rotation matrix from quaternion
C ******************************************************************** C
      ! Returns the 3x3 matrix representing rotation quaternion q
      ! @ input q: Length 4 rotation quaternion
      ! @ returns r: 3x3 rotation matrix
      pure function rot_mat(q) result(r)
      real(8), intent(in) :: q(4)
      real(8)             :: r(3, 3)
      ! Specify the columns
      r(:, 1) = [ one - two * (q(3) ** 2 + q(4) ** 2), 
     1            two * (q(2) * q(3) + q(1) * q(4)), 
     2            two * (q(2) * q(4) - q(1) * q(3)) ]
      r(:, 2) = [ two * (q(2) * q(3) - q(1) * q(4)),
     1            one - two * (q(2) ** 2 + q(4) ** 2),
     2            two * (q(3) * q(4) + q(1) * q(2)) ]
      r(:, 3) = [ two * (q(2) * q(4) + q(1) * q(3)),
     1            two * (q(3) * q(4) - q(1) * q(2)),
     2            one - two * (q(2) ** 2 + q(3) ** 2) ]
      end function
C ******************************************************************** C
C     Stack of skew symmetric matrix from columns of matrix
C ******************************************************************** C
      ! Returns the 9x3 skew symmetric matrix from the columns of R
      ! @ input r: 3x3 matrix
      ! @ outputs rr: Vertical stack of the skew symmetric matrices 
      !               formed from the 3 columns of R
      pure function skew_mat(r) result(rr)
      real(8), intent(in) :: r(3, 3)
      real(8)             :: rr(9, 3)
      !
      rr(1:3, :) = skew_sym(r(:, 1))
      rr(4:6, :) = skew_sym(r(:, 2))
      rr(7:9, :) = skew_sym(r(:, 3))
      end function
C ******************************************************************** C
C     Create a 3x3 matrix out of vector of length 9
C ******************************************************************** C
      ! Returns the 3x3 matrix from the vector v
      ! @ input v: Vector of lenght 9
      ! @ returns r: 3x3 matrix
      pure function vec2mat_9(v) result(r)
      real(8), intent(in) :: v(9)
      real(8)             :: r(3, 3)
      !
      r(:, 1) = [v(1), v(4), v(7)]
      r(:, 2) = [v(2), v(5), v(8)]
      r(:, 3) = [v(3), v(6), v(9)]
      end function
C ******************************************************************** C
C     Compute optimal rotation quaternion
C ******************************************************************** C
      ! Returns the optimal rotation and associated eigenvalue
      ! @ input b: \mathcal{B} matrix from [2]
      ! @ returns lambda_and_q: First entry is the minimum eigenvalue,
      !                         following entries are the associated 
      !                         eigenvector
      function calc_opquat(b) result(lambda_and_q)
      real(8), intent(in) ::  b(4, 4)
      real(8)             ::  be(4, 4)
      real(8)             ::  lambda_and_q(5)
      ! For lapack
      integer             ::  ne, nselect
      parameter               (ne = 4, nselect = 1)
      integer             ::  lda, ldz
      parameter               (lda = ne, ldz = ne)
      integer             ::  info, lwork, liwork, il, iu, m
      real(8)             ::  abstol, vl, vu
      integer :: lwmax
      parameter(lwmax = 1000)
      integer             ::  isuppz(ne), iwork(lwmax)
      real(8)             ::  w(ne), z(ldz, nselect), 
     1                        work(lwmax)
      external dsyevr
      !
      abstol = -1.0
      il = 1
      iu = nselect
      ! lwork and liwork are set based on the query to the optimal 
      ! workspace during testing
      lwork = 104
      liwork = 40
      ! We need B for later, so don't overwrite
      be = b  
      call dsyevr('Vectors', 'Indices', 'Upper', ne, be, lda, vl, vu,il,
     1            iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork,
     2            liwork, info) 
      lambda_and_q(1) = w(1)
      lambda_and_q(2:5) = z(:, 1)
      end function
C ******************************************************************** C
C     Compute G matrix
C ******************************************************************** C
      ! Returns the instantaneous rotation matrix
      ! @ input q: Optimal rotation quaternion
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input b: \matcal{B} matrix from [2]
      ! @ input c: C matrix from [2]
      ! @ input lam: Minimum eigenvalue of B matrix
      ! @ input r: Rotation matrix of optimal rotation quaternion
      ! @ returns g: Instantaneous rotation matrix
      function calc_g(q, n_sh, b, c, lam, r) result(g)
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   ::  q(4), b(4, 4), c(3, n_sh), r(3, 3)
      real(8), intent(in)   ::  lam
      real(8)               ::  g(3, 3*n_sh)
      ! Internal
      real(8)               ::  qrr(4, 4), qrr_3(4, 3), xinv(3, 3), 
     1                          id4(4, 4)
      real(8)               ::  a_temp(3, n_sh)
      ! LAPACK
      integer               :: n_ls, nrhs
      parameter                (n_ls = 3)
      integer               :: lda_ls, ldb_ls
      parameter                (lda_ls = 3, ldb_ls = 3)
      integer               :: info
      integer               :: ipiv(n_ls)
      external dgesv
      ! Function start
      nrhs = 3 * n_sh
      ! Calculate I_4x4
      id4(:, :) = 0.d0
      forall(i = 1:4) id4(i, i) = 1.d0
      ! Calculate instantaneous rotation
      qrr = right_quat(q(1), q(2:4))
      qrr_3 = qrr(:, 2:4)
      xinv = matmul(matmul(transpose(qrr_3), (b - lam * id4)), qrr_3)
      a_temp = matmul(r, c)
      do i = 1, n_sh
        g(:, 3*i-2:3*i) = transpose(-skew_sym(a_temp(:, i)))
      end do
      call dgesv( n_ls, nrhs, xinv, lda_ls, ipiv, g, ldb_ls, info )
      g = 4.d0 * g
      
      end function
C ******************************************************************** C
C     Compute linearized rotation matrix
C ******************************************************************** C
      ! Returns the linearization of the optimal rotation vector
      ! @ input q: Optimal rotation quaternion, size 4
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input g: Instantaneous rotation matrix, 3x3N
      ! @ returns qq: Linearized rotation matrix
      pure function calc_q(q, n_sh, g) result(qq)
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   ::  q(4), g(3, 3*n_sh)
      real(8)               ::  qq(3, 3*n_sh)
      ! Internal
      real(8)               ::  qrr(4, 4), qrr_3(4, 3)
      ! Function start
      ! Compute the linearized rotation matrix
      qrr = right_quat(q(1), q(2:quatdim))
      qrr_3 = qrr(:, 2:quatdim)
      qq = matmul(qrr_3(2:quatdim, :), g)
      
      end function
C ******************************************************************** C
C     Order reference configuration
C ******************************************************************** C
      ! Orders o to have the orientation x, y, z as the columns of o2
      ! @ input o: Un-ordered reference configuration
      ! @ returns o2: Ordered reference configuration such that the 
      !               columns of o2 are a right-handed coordinate system
      pure function order_ini(o) result(o2)
      real(8), intent(in) :: o(3, 3)
      real(8)             :: o2(3, 3), xo(3), yo(3), zo(3), cross_vec(3)
      real(8)             :: test
      ! z is the min eigenvalue, y is the greatest eigenvalue, 
      ! x is the remaining eigenvalue
      zo = o(:, 1)
      yo = o(:, 3)
      xo = o(:, 2)
      ! Test to see if right-handed system
      cross_vec(1) = yo(2) * zo(3) - yo(3) * zo(2)
      cross_vec(2) = -(yo(1) * zo(3) - yo(3) * zo(1))
      cross_vec(3) = yo(1) * zo(2) - yo(2) * zo(1)
      test = dot_product(xo, cross_vec)
      if (test > 0) then
        ! Correct assumption, right hand system
        o2(:, 1) = xo
        o2(:, 2) = yo
        o2(:, 3) = zo
      else
        ! Reverse the sense of z to make it right-handed
        o2(:, 1) = xo
        o2(:, 2) = yo
        o2(:, 3) = -zo
      end if
      end function
C ******************************************************************** C
C     Compute reference configuration
C ******************************************************************** C
      ! Returns the right-hand ordered reference configuration
      ! @ input init_pts: initial [x, y, z] coordinates of each point
      ! @ input n_sh: Number of shell nodes (N)
      ! @ returns o: Orientation of the reference configuration that is 
      !              a valid right-handed coordinate system.
      function ini_config(init_pts, n_sh) result(o)
      ! Input and output
      integer, intent(in) ::  n_sh
      real(8), intent(in) ::  init_pts(3, n_sh)
      real(8)             ::  o(3, 3)
      ! Internal
      integer             :: ne, nselect
      parameter             (ne = 3, nselect = 3)
      integer             :: lda, ldz
      parameter             (lda = ne, ldz = ne)
      integer             :: info, lwork, liwork, il, iu, m
      real(8)             :: abstol, vl, vu
      integer             :: lwmax
      parameter             (lwmax = 1000)
      integer             :: isuppz(ne), iwork(lwmax)
      double precision    :: ae(lda, ne), w(ne), z(ldz, nselect), 
     1                       work(lwmax)
      external dsyevr
      !
      ae(:, :) = zero
      ! The spread is used to compute the outer product of two vectors
      do i = 1, n_sh
        ae = ae 
     1  + spread(init_pts(:, i), dim=2, ncopies=3) *
     2    spread(init_pts(:, i), dim=1, ncopies=3)
      end do
      !
      abstol = -1.0
      il = 1
      iu = nselect
      ! lwork and liwork are set based on the query to the optimal 
      ! workspace during testing
      lwork = 104
      liwork = 40
      call dsyevr('Vectors', 'A', 'Upper', ne, ae, lda, vl, vu, il,
     1            iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork,
     2            liwork, info) 
      ! The eigenvalues are returned in asscending order
      o = order_ini(z)
      end function
C ******************************************************************** C
C     Calculate the warping function for all nodes
C ******************************************************************** C
      ! Returns the warping function evaluated at each node.
      ! @ input c: Location of nodes in initial config. relative to 
      !            centroid, size 3xN
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input r0: Rotation from reference to initial configuration
      ! @ returns psi: Warping function at each node.
      !
      ! Notes:
      !   - If the distance between adjacent nodes is less than 1e-8
      !   then nodes on the web can be incorrectly identified
      pure function warp_fun(c, n_sh, r0) result(psi)
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   :: c(3, n_sh), r0(3, 3)
      real(8)               :: psi(n_sh), c2(3, n_sh)
      ! Internal
      real(8)             :: s, s1, s2, two2, four, tol
      parameter              (two2 = 2.d0, four = 4.d0, tol = 1.d-8)
      integer             :: ii
      ! Allocate array
      ! Rotate from the initial config to the reference config
      c2 = matmul(transpose(r0), c)
      ! Calculate the width and depth of the section
      s1 = two2 * maxval(c2(1, :))
      s2 = two2 * maxval(c2(2, :))
      ! Calculate the warping function, \psi = X*Y
      ! todo: clean this function up
      do ii = 1, n_sh
        if (abs(c2(2, ii) - s2 / two2 ) < tol) then
          ! Node is on the top flange
          s = c2(1, ii) + s1 / two2
          !psi(ii) = -s2 * (s1 - two2 * s) / four
          psi(ii) = c2(1, ii) * c2(2, ii)
        else if (abs(c2(2, ii) + s2 / two2 ) < tol) then
          ! Node is on the bottom flange
          s = c2(1, ii) + s1 / two2
          !psi(ii) = s2 * (s1 - two2 * s) / four
          psi(ii) = c2(1, ii) * c2(2, ii)
        else
          ! Node is on the web line
          psi(ii) = 0.d0
        end if
      end do
      end function
C ******************************************************************** C
C     Calculate the warping amplitude
C ******************************************************************** C
      ! Returns the value of the warping amplitude.
      ! @ input x_ini: pos. of each node in the initial configuration
      ! @ input x_def: pos. of each node in the deformed configuration
      ! @ input u_c: Translation of the centroid
      ! @ input t: Orientation of the normal to the cross-section
      ! @ input r: Rotation from initial to deformed configurations
      ! @ input psi: Warping function at each node
      ! @ returns w: Warping amplitude
      pure function warp_amp(x_ini, x_def, u_c, t, r, psi) result(w)
      ! Input and output
      real(8), intent(in) ::  x_def(:, :), x_ini(:, :), u_c(3), t(3)
      real(8), intent(in) ::  r(3, 3), psi(:)
      real(8)             ::  w, resid(3)
      ! Internal
      integer             :: ii, non_zero_nodes, sz(2)
      real(8)             :: zero_tol
      parameter              (zero_tol=1.d-6)
      ! Function start
      sz = shape(x_def)
      w = 0.d0
      non_zero_nodes = 0
      do ii = 1, sz(2)
        if (abs(psi(ii)) > zero_tol) then
          non_zero_nodes = non_zero_nodes + 1
          resid = x_def(1:3, ii) - matmul(r, x_ini(1:3, ii) + u_c)
          w = w + dot_product(t, resid) / psi(ii)
        end if
      end do
      w = w / non_zero_nodes
      end function
C ******************************************************************** C
C     Compute the linearized warping vector
C ******************************************************************** C
      ! todo: decide what form we want and clean this function up
      ! Returns the linearized warping vector.
      ! @ input x_def: pos. of nodes in deformed configuration, 3xN
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input t0: Orientation of the normal to the cross-section in  
      !             the reference configuration
      ! @ input psi: Warping function for all the nodes
      ! @ input r: Rotation from initial to deformed configuration
      ! @ input g: Instantaneous rotation matrix, size 3x3N
      ! @ input we: Area weight value for each node, size N
      ! @ returns w_lin: Linearization of warping amplitude w.r.t x
      !
      ! Notes:
      !   - we containts the area weights for each node, then the shear
      !   weights for each node, only use the first set
      pure function calc_lin_w(c2, n_sh, t0, psi, r, we) 
     1              result(w_lin)
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   ::  c2(3, n_sh), t0(3), psi(n_sh),r(3, 3)
      real(8), intent(in)   ::  we(n_sh)
      real(8)               ::  w_lin(3*n_sh)
      ! Internal
      integer               ::  i, j, num_psi_non_zero
      real(8)               ::  t(3), we_total
      real(8)               ::  zero, one, two
      parameter                 (zero = 0.d0, one = 1.d0, two = 2.d0)
      real(8)               ::  zero_tol, d_cl
      parameter                 (zero_tol=1.d-6)
      ! Assemble the linearized warping vector
      t = matmul(r, t0)
      w_lin = 0.d0
      we_total = 0.d0
      do i = 1, n_sh
        if (abs(psi(i)) > zero_tol) then
          ! todo: bimom can be calculated in the function with psi
          we_total = we_total + abs(c2(1, i) * psi(i) * we(i))
          w_lin(3*i-2:3*i) = w_lin(3*i-2:3*i) + t * psi(i) * we(i)
        end if
      end do
      ! Normalize by the total bimoment to satisfy equilibrium
      ! B = (M_f,t + M_f,b) * (d_cl / 2)
      d_cl = maxval(c2(2, :))
      we_total = we_total * d_cl
      w_lin = w_lin / we_total      
      
      end function
C ******************************************************************** C
C     Extract rotation vector from quaternion
C ******************************************************************** C
      ! Returns the rotation vector extracted from the quaternion
      ! @ input q: Rotation quaternion
      ! @ returns phi: Euler angle representation of q
      !
      ! From Abaqus Theory Guide 1.3.1
      pure function extract_rotation(q) result(phi)
      ! Input and output
      real(8), intent(in) ::  q(4)
      real(8)             ::  phi(3)
      ! Internal
      real(8)             :: q_mag, rot
      real(8)             :: small
      parameter              (small = 1.d-14)
      ! Function start
      q_mag = sqrt(q(2) ** 2 + q(3) ** 2 + q(4) ** 2)
      if (q_mag > small) then
        rot = 2.d0 * atan2(q_mag, q(1))
        phi = rot / q_mag * q(2:4)
      else
        phi(:) = 0.d0
      end if
      end function
C ******************************************************************** C
C     Determine the weights for all the nodes
C ******************************************************************** C
      ! Returns the weight value for each node
      ! @ input p: Shell nodes in the reference configuration, size 3xN
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input tf_tw_code: Encoded flange and web thickness factors
      ! @ input etang: Tangent modulus, size N
      ! @ returns w: Weight value for each node
      !   - Row 1: Area weight
      !   - Row 2: Normalized strong axis weight
      !   - Row 3: Normalized weak axis weight
      !   - Row 4: Normalized complementary strong axis weight
      !
      ! Notes:
      !   - The weight represents the equivalent area of each node
      !   - tf_tw_code is defined as [d][tf][tw], where d is the number
      !     of decimals in tf and tw
      !   - The weights are only valid if the interface is elastic
      !   - The flange and web mesh distances are assumed to be constant
      !     between nodes
      function calc_weights(p, n_sh, tf_tw_code, etang) result(ww)
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   ::  p(3, n_sh), etang(n_sh)
      integer, intent(in)   ::  tf_tw_code
      real(8)               ::  ww(4, n_sh)
      ! Internal variables
      ! Tags for node positions
      integer               ::  web_tag, flange_tag, corner_tag, 
     1                          joint_tag
      parameter                 (web_tag=1, flange_tag=2, corner_tag=3,
     1                           joint_tag=4)
      integer               ::  classification(n_sh)
      ! For weights at each tag
      real(8)               ::  web_weight, flange_weight,corner_weight,
     1                          joint_weight
      real(8)               ::  tf, tw, delta_f, delta_w, x_max, y_max,
     1                          pt(3), x_max_2, y_max_2
      integer               ::  ii, d_digits, deci, n_digits, c, tf_tw
      ! For shear weights
      real(8) ::  width, depth, v_weights(n_sh, 3)
      ! Tolerance in positions
      real(8) :: zero_tol
      parameter(zero_tol=1.d-3)
      
      ! Function start
      ! Classify the points
      classification(:) = 0
      x_max = maxval(p(1, :))
      y_max = maxval(p(2, :))
      do ii = 1, n_sh
        pt = p(:, ii)
        if (abs(pt(1)) < zero_tol) then
          ! Point is on the web
          classification(ii) = web_tag
        end if
        if (abs(abs(pt(2)) - y_max) < zero_tol) then
          ! Point is on the flange
          if (classification(ii) == web_tag) then
            ! Point is on the web and flange
            classification(ii) = joint_tag;
          else if (abs(abs(pt(1)) - x_max) < zero_tol) then
            ! Point is on the corner of the flange
            classification(ii) = corner_tag
          else
            ! Just on the flange
            classification(ii) = flange_tag
          end if
        end if
      end do
      
      ! Calculate the mesh distance for flange and web nodes
      ! Start search for 2nd greatest values
      x_max_2 = -x_max
      y_max_2 = -y_max
      do ii = 1, n_sh
        pt = p(:, ii)
        c = classification(ii)
        ! Web
        if (c == web_tag) then
          ! Set 2nd largest y value
          if (y_max_2 < pt(2) .and. pt(2) < y_max) then
            y_max_2 = pt(2)
          end if
        end if
        ! Flange
        if ((c==flange_tag .or. c==joint_tag) .and. pt(2) > 0.d0) then
          ! Set the 2nd largest x value
          if (x_max_2 < pt(1) .and. pt(1) < x_max) then
            x_max_2 = pt(1)
          end if
        end if
      end do
      delta_f = x_max - x_max_2
      delta_w = y_max - y_max_2
      
      ! Decode the flange and web thickness factors
      ! d_digits starts at 4 since there will be at least 5 digits
      d_digits = 4;
      do while (tf_tw_code / 10 ** d_digits > 10)
          ! Divide by 10 until we get something less than 10, this is 
          ! the first digit of tf_tw_code
          d_digits = d_digits + 1
      end do
      deci = tf_tw_code / 10 ** d_digits
      ! The number of digits in each thickness is assumed to be equal
      n_digits = d_digits / 2
      ! The thickness values are what is left after we subtract the 
      ! decimal value from the front
      tf_tw = tf_tw_code - deci * (10 ** d_digits)
      ! Use the same subtracting off-the-front to get the tf and then tw
      ! The int casting is to truncate the un-needed digits
      tf = int(tf_tw / 10 ** n_digits)
      tf = tf / 10. ** deci
      tw = int(tf_tw - int(tf * 10. ** (deci + n_digits)))
      tw = tw / 10. ** deci
      
      ! Calculate the area weights for each class
      web_weight = delta_w * tw
      flange_weight = delta_f * tf
      corner_weight = 0.5d0 * flange_weight
      joint_weight = flange_weight + 0.5d0 * tw * delta_w
      
      ! Calculate the shear weights
      width = 2.d0 * x_max
      depth = 2.d0 * y_max + tf
      v_weights = shear_factors(p, n_sh, depth, width, tf, tw, 
     1                          delta_f, delta_w, classification, etang)
      ! Assign the area weights
      do ii = 1, n_sh
        c = classification(ii)
        if (c == web_tag) then
          ww(1, ii) = web_weight
        else if (c == flange_tag) then
          ww(1, ii) = flange_weight
        else if (c == corner_tag) then
          ww(1, ii) = corner_weight
        else if (c == joint_tag) then
          ww(1, ii) = joint_weight
        end if
        ! Add the shear weights
        ww(2:4, ii) = v_weights(ii, :)
      end do
      end function
C ******************************************************************** C
C     Calculate the shear factor for flange nodes (weak and strong)
C ******************************************************************** C
      ! Returns the shear factors for the flange nodes
      ! @ input xi: Flange node
      ! @ input d: Section depth
      ! @ input d1: Section depth minus two flange thicknesses
      ! @ input bf: Section flange width
      ! @ input tf: Flange thickness
      ! @ input tw: Web thickness
      ! @ returns: The shear resultants
      !             v(1): strong axis flange
      !             v(2): weak axis flange
      pure function shear_flange(xi, d, d1, bf, tf, tw) result(v)
        ! Input and output
        real(8), intent(in) ::  xi(3), d, d1, bf, tf, tw
        real(8)             ::  v(2)
        ! Internal
        real(8)             ::  x1, y1
        ! Function start
        x1 = xi(1)
        y1 = xi(2)
        ! Strong axis
        v(1) = (d ** 2 - 4. * y1 ** 2) / 8.
        ! Weak axis
        v(2) = (bf ** 2 - 4. * x1 ** 2) / 8.
      end function
C ******************************************************************** C
C     Calculate the complementary shear factor for flange nodes 
C ******************************************************************** C
      ! Returns the complementary shear factors for the flange nodes
      ! @ input xi: Flange node
      ! @ input d_cl: Section depth minus flange thicknesses
      ! @ input bf: Section flange width
      ! @ returns: The complementary shear in the direction of the weak
      !            axis
      !
      ! Notes:
      !   - This function assumes that the force is applied in the posi-
      !   tive y direction
      !   - At the joint nodes no complementary shear is applied if the
      !   positive and negative forces cancel out
      pure function comp_shear_flange(xi, d_cl, bf) result(v)
        ! Input and output
        real(8), intent(in) ::  xi(3), d_cl, bf
        real(8)             ::  v
        ! Internal
        real(8)             ::  x1, y1, sgn
        ! Function start
        x1 = xi(1)
        y1 = xi(2)
        ! Direction assumes a positive shear force
        sgn = sign(1.d0, x1 * y1)
        ! Linear interpolation based on position along half flange
        if (abs(xi(1)) > bf / 2.) then
          v = 0.
        else
          v = sgn * (bf - 2. * abs(x1)) * d_cl / 4.
        end if
      end function
C ******************************************************************** C
C     Calculate the shear factor for web nodes (weak and strong)
C ******************************************************************** C
      ! Returns the shear factors for the web nodes
      ! @ input xi: Flange node
      ! @ input d: Section depth
      ! @ input d1: Section depth minus two flange thicknesses
      ! @ input bf: Section flange width
      ! @ input tf: Flange thickness
      ! @ input tw: Web thickness
      ! @ returns: The shear resultants
      !             v(1): strong axis web 
      !             v(2): weak axis web 
      pure function shear_web(xi, d, d1, bf, tf, tw) result(v)
        ! Input and output
        real(8), intent(in) ::  xi(3), d, d1, bf, tf, tw
        real(8)             ::  v(2)
        ! Internal
        real(8)             ::  x1, y1
        ! Function start
        x1 = xi(1)
        y1 = xi(2)
        ! Strong axis
        if (abs(xi(2)) > (d - tf) / 2.) then
          v(1) = 0.
        else 
          v(1) = (2.*bf*(d+d1)*tf + tw*(d1**2-4.*y1**2)) / (8.*tw)
        end if
        ! Weak axis
        ! Only valid for x1 = 0
        v(2) = (2.*bf**2*tf + (d-2.*tf)*tw**2) / (8.*d1)
      end function
C ******************************************************************** C
C     Calculate the shear factor resultants for the strong and weak axes
C ******************************************************************** C
      ! Returns the shear factor resultant for one flange and web
      ! @ input p: Shell nodes in the reference configuration, 3xN
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input d: Section depth
      ! @ input bf: Section width
      ! @ input tf: Flange thickness
      ! @ input tw: Web thickness
      ! @ input c_tags: Node classification tag, size N
      ! @ input etang: Tangent modulus, size N
      ! @ returns: The normalized shear resultants, size Nx3
      !   - Column 1: strong axis
      !   - Column 2: weak axis
      !   - Column 3: complementary force for strong axis
      pure function shear_factors(p, n_sh, d, bf, tf, tw, 
     1                            delta_f, delta_w, c_tags, etang)
     2                            result(v)
      ! Input and output
      integer, intent(in) ::  n_sh, c_tags(n_sh)
      real(8), intent(in) ::  p(3, n_sh), d, bf, tf, tw, 
     1                        delta_f, delta_w, etang(n_sh)
      real(8)             ::  v(n_sh, 3)
      ! Internal variables
      real(8)             ::  d1, tau(2), total_strong, total_weak,
     1                        p_corner(3), tau_2, d_cl
      integer             ::  ii, c
      ! Classification tags
      integer             ::  web_tag, flange_tag, corner_tag, 
     1                        joint_tag
      parameter              (web_tag=1, flange_tag=2, corner_tag=3,
     1                        joint_tag=4)
      ! Function start
      ! The first value of tau is for the strong axis, the second is for
      ! the weak axis
      ! tau_2 is the complementary shear in the flange
      d1 = d - 2. * tf
      d_cl = d - tf
      total_strong = 0.
      total_weak = 0.
      do ii = 1, n_sh
        c = c_tags(ii)
        if (c == web_tag) then
          tau = tw * delta_w * shear_web(p(:, ii), d, d1, bf, tf, tw)
          ! Don't consider the complementary shear for the web
          tau_2 = 0.
        else if (c == flange_tag) then
          tau = tf * delta_f * shear_flange(p(:, ii), d, d1, bf, tf, tw)
          tau_2 = tf * delta_f * comp_shear_flange(p(:, ii), d_cl, bf)
          
          tau(1) = 0.d0
        else if (c == corner_tag) then
          ! Sample at one-half element distance towards the web
          p_corner(1) = -sign(1., p(1, ii)) * delta_f / 2.
          p_corner(2:3) = 0.
          p_corner = p(:, ii) + p_corner
          tau = tf*delta_f/2.* shear_flange(p_corner, d, d1, bf, tf, tw)
          tau_2 = tf * delta_f / 2.* comp_shear_flange(p_corner,d_cl,bf)
          
          tau(1) = 0.d0
        else
          tau=tf*delta_f*shear_flange(p(:, ii),d,d1,bf,tf,tw)
          
          tau(1) = 0.d0
          
          tau = tau+tw*delta_w/2.*shear_web(p(:, ii), d, d1, bf, tf, tw)
          ! Complementary shear cancels from both sides at the joint
          tau_2 = 0.
        end if
        ! Adjust for the tangent modulus
        tau = tau * etang(ii)
        ! Record the weights
        total_strong = total_strong + tau(1)
        total_weak = total_weak + tau(2)
        v(ii, 1:2) = tau
        v(ii, 3) = tau_2
      end do
      
      ! Normalize the weights, the complementary takes the same 
      ! normalization as the strong axis
      v(:, 1) = v(:, 1) / total_strong
      v(:, 2) = v(:, 2) / total_weak
      v(:, 3) = v(:, 3) / total_strong
      end function
C ******************************************************************** C
C     Calculate the linearized displacement matrix
C ******************************************************************** C
      ! Returns the matrix relating centroidal disp. to shell disps.
      ! @ input w_all: Area and shear weights, size 4xN
      ! @ input a_total: Total area
      ! @ input r0: Initial orientation, size 3x3
      ! @ input r: Rotation from initial to deformed, size 3x3
      ! @ returns: The displacement factor for each direction for each 
      !            node, size 3xN
      pure function calc_disp_lin(n_sh, w_all, a_total, r0, r) 
     1              result(u_lin)
      ! Input and output
      integer, intent(in) ::  n_sh
      real(8), intent(in) ::  w_all(4, n_sh), a_total, r0(3, 3), r(3, 3)
      real(8)             ::  u_lin(3, 3 * n_sh)
      ! Internal variables
      real(8)             ::  w_bar(3, 3), rr0(3, 3), s_bar(3, 3), 
     1                        t_bar(3), s(3, 3), rr0_T(3, 3)
      integer             ::  ii   
      
      ! Function start
      t_bar = r0(:, 3)
      rr0 = matmul(r, r0)
      rr0_T = transpose(rr0)
      do ii = 1, n_sh
        ! Weak axis force due to weak axis displacement
        w_bar = 0.d0
        w_bar(1, 3) = w_all(3, ii)
        w_bar(3, 1) = w_all(3, ii)
        s_bar(:, 1) = matmul(w_bar, t_bar)
        ! Weak axis force in flange due to strong axis displacement
        w_bar = 0.d0
        w_bar(3, 1) = w_all(4, ii)
        w_bar(1, 3) = w_all(4, ii)
        ! Strong axis force due to strong axis displacement
        w_bar(3, 2) = w_all(2, ii)
        w_bar(2, 3) = w_all(2, ii)
        s_bar(:, 2) = matmul(w_bar, t_bar)
        ! Axial force due to axial displacement
        w_bar = 0.d0
        w_bar(3, 3) = w_all(1, ii) / a_total
        s_bar(:, 3) = matmul(w_bar, t_bar)
        ! Rotate the weights from the reference to the deformed config
        s = matmul(rr0, matmul(s_bar, rr0_T))
        u_lin(:, 3*ii-2:3*ii) = transpose(s)
      end do
      end function
C ******************************************************************** C
C     Correct tangential stress under bending
C ******************************************************************** C
      ! Returns the matrix of tangential forces for bending axial stress
      ! @ input 
      ! @ returns: 
      pure function tang_force_b(n_sh, p, r0, r, df, dw, tf, tw,b,hcl) 
     1              result(f_tang_b)
      ! Input and output
      integer, intent(in) ::  n_sh
      real(8), intent(in) ::  p(3, n_sh), r0(3, 3), r(3, 3), df, dw, tf, 
     1                        tw, b, hcl
      real(8)             ::  f_tang_b(3, n_sh)
      ! Internal variables
      integer             ::  i
      real(8)             ::  moi,nu,x_vec(3),y_vec(3),sf,zero_tol,y
      parameter               (sf = 1.0d0, zero_tol = 0.1)
      
      ! Assume Poisson ratio = 0.3
      nu = 0.3
      ! moi is the moment of intertia
      moi = 0.d0
      x_vec = matmul(r, r0(:, 1))
      y_vec = matmul(r, r0(:, 2))
      f_tang_b = 0.d0
      do i = 1, n_sh
        ! web
        y = p(2, i)
        if (abs(p(1, i)) < zero_tol) then
          if (abs(abs(p(2, i)) - hcl / 2.) < zero_tol) then
            ! joint nodes
         f_tang_b(2, i)=web_f_r(y+dw/2.,dw,hcl)-web_f_r(y-dw/2.,dw,hcl)
            f_tang_b(2, i) = f_tang_b(2, i) * nu * tw * dw / 2.d0
            moi = moi + p(2, i) ** 2 * tw * dw / 2.
          else
            ! inner nodes
            if (abs(p(2, i)) < 0.1) then
            ! center node (now nothing different from other nodes)
         f_tang_b(2, i)=web_f_r(y+dw/2.,dw,hcl)-web_f_r(y-dw/2.,dw,hcl)
            f_tang_b(2, i) = f_tang_b(2, i) * nu * tw * dw / 2.d0
            else
         f_tang_b(2, i)=web_f_r(y+dw/2.,dw,hcl)-web_f_r(y-dw/2.,dw,hcl)
            f_tang_b(2, i) = f_tang_b(2, i) * nu * tw * dw / 2.d0
            end if
            moi = moi + p(2, i) ** 2 * tw * dw
          end if
          ! sf = scale factor to account for flexibility
          f_tang_b(:, i) = f_tang_b(2, i) * y_vec * sf
        end if
        ! flange
        if (abs(abs(p(2, i)) - hcl / 2.) < zero_tol) then
          if (abs(abs(p(1, i)) - b / 2.) < zero_tol) then
            ! edge node
            moi = moi + p(2, i) ** 2 * tf * df / 2.
          else
            ! inner node
            moi = moi + p(2, i) ** 2 * tf * df
          end if
        end if
      end do
      
      ! Normalize to unit moment
      f_tang_b = f_tang_b / moi
      end function      
C ******************************************************************** C
      ! Returns restraint force
      ! @ input 
      ! @ returns: 
      pure function web_f_r(y, dw, hcl) result(fdiff)
      ! Input and output
      real(8), intent(in) ::  y, dw, hcl
      real(8)             ::  fdiff
      real(8)             ::  alpha, beta
      ! Function
      alpha = 0.10
      beta = 2.0
      if (abs(y) <= hcl / 2.d0) then
        fdiff = -y / (alpha * abs(y) / (hcl / 2.d0) + 1.d0 + beta)
      else
        fdiff = 0.d0
      end if
      end function
C ******************************************************************** C
      ! Returns torque on flange only
      ! @ input 
      ! @ returns: 
      function torquer(n_sh, p_ref, rr, we, hcl) result(tq)
      ! Input and output
      integer, intent(in) ::  n_sh
      real(8), intent(in) ::  p_ref(3, n_sh), rr(3, 3), we(n_sh), hcl
      real(8)             ::  tq(3, n_sh)
      real(8)             ::  tor, ff
      integer             ::  ii
      ! Function
      
      ! todo: needs the initial config, r0, then apply r . r0 . tq
      
      tor = 0.
      tq = 0.d0
      do ii = 1, n_sh
        if (abs(p_ref(2, ii)) >= hcl / 2.d0) then
          ! todo: fix this hack (the 150 and 0.5 factor) - concentrates force at the center nodes
          ff = 1.d0 - 0.9 * abs(p_ref(1, ii)) / 150.
          tq(1:2, ii) = ff * we(ii) * [ -p_ref(2, ii), p_ref(1, ii) ]
          tq(1:3, ii) = matmul(rr, tq(1:3, ii))
          tor = tor + ff * we(ii) * norm2(p_ref(1:2, ii)) ** 2
        end if
      end do
      tq = -1.d0 * tq / tor
      end function
C ******************************************************************** C
      pure function rot_weights(n_sh, xyz, area, tmod) result(rw)
        ! Returns the weights for applied moments.
        ! @input n_sh: Number of shell nodes on the interface
        ! @input xyz: Coords of shell nodes in reference config
        ! @input area: Tributary area of each shell node
        ! @input tmod: Tangent modulus associated with each node
        ! @returns: Mx,My,Mz weights
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  xyz(3, n_sh), area(n_sh), tmod(n_sh)
        integer             ::  ii
        real(8)             ::  rw(4, n_sh)
        real(8)             ::  mx(n_sh), my(n_sh), mz(2, n_sh), mxtot,
     1                          mytot, mztot
        do ii = 1, n_sh
          mx(ii) = xyz(2, ii) * area(ii) * tmod(ii)
          mxtot = mxtot + mx(ii) * xyz(2, ii)
          my(ii) = -xyz(1, ii) * area(ii) * tmod(ii)
          mytot = mytot + my(ii) * xyz(1, ii)
          mz(1:2, ii) = [-xyz(2, ii), xyz(1, ii)] * area(ii) * tmod(ii)
          mztot = mztot + area(ii) * tmod(ii) * norm2(xyz(1:2, ii))**2
        end do
        rw(1, :) = mx / mxtot
        rw(2, :) = my / abs(mytot)
        rw(3:4, :) = mz / mztot
      end
        
C ******************************************************************** C
      pure function lin_rot(n_sh, r0, r, rw) result(linr)
        ! Returns the linearization of the rotation constraints.
        ! @input n_sh: Number of shell nodes on the interface
        ! @input r0: Rotation matrix from reference to initial configs
        ! @input r: Rotation matrix from initial to current configs
        ! @input rw: Weights for the moments
        ! @returns: Matrix for the linearization of the rotations
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  r(3, 3), r0(3, 3), rw(4, n_sh)
        real(8)             ::  linr(3, 3*n_sh)
        integer             ::  ii
        real(8)             ::  w_bar(3, 3), s_bar(3, 3), s(3, 3),
     1                          t_bar(3), rr0(3, 3), rr0_T(3, 3)
        t_bar = r0(:, 3)
        rr0 = matmul(r, r0)
        rr0_T = transpose(rr0)
        do ii = 1, n_sh
          ! Moment about x
          w_bar = 0.d0
          w_bar(3, 3) = rw(1, ii)
          s_bar(:, 1) = matmul(w_bar, t_bar)
          ! Moment about y
          w_bar = 0.d0
          w_bar(3, 3) = rw(2, ii)
          s_bar(:, 2) = matmul(w_bar, t_bar)
          ! Moment about z (torque)
          w_bar = 0.d0
          w_bar(1, 3) = rw(3, ii)
          w_bar(3, 1) = rw(3, ii)
          w_bar(2, 3) = rw(4, ii)
          w_bar(3, 2) = rw(4, ii)
          s_bar(:, 3) = matmul(w_bar, t_bar)
          ! Rotate from reference to current config
          s = matmul(rr0, matmul(s_bar, rr0_T))
          linr(1:3, 3*ii-2:3*ii) = transpose(s)
        end do        
      end
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
      end  ! END SUBROUTINE
