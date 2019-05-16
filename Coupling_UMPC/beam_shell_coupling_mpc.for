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
      real(8)               :: total_weight
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
     3                         w_lin(:), psi_all(:), weights(:, :)
      ! For the rotation quaternion and computations
      real(8)               :: b_mat(quatdim, quatdim), 
     1                         beta(quatdim, quatdim)
      real(8)               :: op_quat(quatdim)
      real(8)               :: lambda
      real(8)               :: lam_q(quatdim+1)
      ! Linearized displacment
      real(8), allocatable  ::  disp_lin(:, :)
      ! Numerical parameters
      real(8)               :: one, two, zero
      parameter                (one=1.0d0, two=2.0d0, zero=0.d0)
      ! Allocate variable size arrays
      n_shell = n - 1
      allocate(x_shell(ndim, n_shell))
      allocate(u_shell(ndim, n_shell))
      allocate(x_shell_def(ndim, n_shell))
      allocate(c_mat(ndim, n_shell))
      allocate(d_mat(ndim, n_shell))
      allocate(g_mat(ndim, ndim * n_shell))
      allocate(q_mat(ndim, ndim * n_shell))
      allocate(psi_all(n_shell))
      allocate(w_lin(ndim * n_shell))
      allocate(weights(4, n_shell))
      allocate(disp_lin(ndim, ndim * n_shell))
C ******************************************************************** C
C Subroutine definition
C ******************************************************************** C
      ! Calculate the centroid locations
      ! The first node is the beam node, the shell nodes follow
      x_shell = x(1:3, 2:n)
      u_shell = u(1:3, 2:n)
      x_shell_def = x_shell + u_shell
      xc0(:) = zero
      xc(:) = zero
      do i = 1, n_shell
        xc0(:) = xc0(:) + x_shell(:, i)
        xc(:) = xc(:) + x_shell_def(:, i)
      end do
      n_reciprocal = one / n_shell
      xc0 = xc0 * n_reciprocal
      xc = xc * n_reciprocal
      u_cent = xc - xc0
      
      ! Find the reference orientation
      do i = 1, n_shell
        c_mat(:, i) = x_shell(:, i) - xc0
        d_mat(:, i) = x_shell_def(:, i) - xc
      end do
      r_ref = ini_config(c_mat, n_shell)
      
      ! Calculate the nodal weights
      ! The centroid does not depend on the weights due to symmetry so
      ! it's valid to calculate the centroid first.
      ! The first n_shell values are the area weights, the rest are 
      ! shear weights.
      weights = calc_weights(matmul(transpose(r_ref), c_mat), n_shell, 
     1                       jtype)
      total_weight = zero
      do i = 1, n_shell
        total_weight = total_weight + weights(1, i)
      end do
      
      ! Compute the optimal rotation quaternion
      b_mat(:, :) = zero
      do i = 1, n_shell
        beta = left_quat(zero,d_mat(:, i))- right_quat(zero,c_mat(:, i))
        b_mat = b_mat + weights(1, i) * matmul(transpose(beta), beta)
      end do
      lam_q = calc_opquat(b_mat)
      lambda = lam_q(1)
      op_quat = lam_q(2:5)
      ! Make the first entry always positive
      op_quat = sign(one, op_quat(1)) * op_quat  
      
      ! Compute the linearized rotation matrix
      r_mat = rot_mat(op_quat)
      g_mat = calc_g(op_quat, n_shell, b_mat, c_mat, lambda, r_mat)
      q_mat = calc_q(op_quat, n_shell, g_mat)
      
      ! Compute the warping amplitude
      psi_all = warp_fun(c_mat, n_shell, r_ref)
      t_def = matmul(r_mat, r_ref(:, 3))
      w_amp = warp_amp(x_shell, x_shell_def, u_cent, t_def, r_mat,
     1                 psi_all)
     
      ! Compute the linearized warping
      w_lin = calc_lin_w(c_mat, n_shell, psi_all, r_ref, r_mat, 
     1                   weights(1, :))
      
      ! Compute the linearized centroid displacement
      disp_lin = calc_disp_lin(n_shell, weights, total_weight, r_ref, 
     1                         r_mat)
      ! Assign the A submatrix for the beam node and set active DOF
      forall(i = 1:maxdof) a(i, i, 1) = one
      forall(i = 1:maxdof) jdof(i, 1) = i
      
      ! Assign the A submatrices for the shell nodes
      do i = 2, n
        ! Index corresponding to the shell nodes
        i_sh = i - 1
        ! Displacement constraints
        a(1:3, 1:3, i) = -disp_lin(1:3, 3*i_sh-2:3*i_sh)
        ! Rotation constraints
        a(4, 1:3, i) = -q_mat(1, 3*i_sh-2:3*i_sh) * weights(1, i_sh)
        a(5, 1:3, i) = -q_mat(2, 3*i_sh-2:3*i_sh) * weights(1, i_sh)
        a(6, 1:3, i) = -q_mat(3, 3*i_sh-2:3*i_sh) * weights(1, i_sh)
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
      integer             :: i
      ! 
      ql(:, 1) = [ v0, v(:) ]
      ql(1, 2:4) = -v(:)
      ql(2:4, 2:4) = skew_sym(v)
      forall(i = 1:3) ql(1 + i, 1 + i) = ql(1 + i, 1 + i) + v0
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
      integer             :: i
      !
      qr(:, 1) = [ v0, v(:) ]
      qr(1, 2:4) = -v(:)
      qr(2:4, 2:4) = -skew_sym(v)   
      forall(i = 1:3) qr(1 + i, 1 + i) = qr(1 + i, 1 + i) + v0   
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
      real(8)             ::  ae(lda, ne), w(ne), z(ldz, nselect), 
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
      integer               ::  sz(2)
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
      integer               ::  sz(2)
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
      integer             :: ne, nselect, mat_shape(2)
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
      real(8)             :: s, s1, s2, two, four, tol
      parameter              (two = 2.d0, four = 4.d0, tol = 1.d-8)
      integer             :: i
      ! Allocate array
      ! Rotate from the initial config to the reference config
      c2 = matmul(transpose(r0), c)
      ! Calculate the width and depth of the section
      s1 = two * maxval(c2(1, :))
      s2 = two * maxval(c2(2, :))
      ! Calculate the warping function, \psi = X*Y
      do i = 1, n_sh
        if (abs(c2(2, i) - s2 / two ) < tol) then
          ! Node is on the top flange
          s = c2(1, i) + s1 / two
          psi(i) = -s2 * (s1 - two * s) / four
        else if (abs(c2(2, i) + s2 / two ) < tol) then
          ! Node is on the bottom flange
          s = c2(1, i) + s1 / two
          psi(i) = s2 * (s1 - two * s) / four
        else
          ! Node is on the web line
          psi(i) = 0.d0
        end if
      end do
      end function
C ******************************************************************** C
C     Calculate the warping amplitude
C ******************************************************************** C
      ! Returns the value of the warping amplitude.
      ! @ input x_ref: pos. of each node in the initial configuration
      ! @ input x_def: pos. of each node in the deformed configuration
      ! @ input u_c: Translation of the centroid
      ! @ input t: Orientation of the normal to the cross-section
      ! @ input r: Rotation from initial to deformed configurations
      ! @ input psi: Warping function at each node
      ! @ returns w: Warping amplitude
      pure function warp_amp(x_ref, x_def, u_c, t, r, psi) result(w)
      ! Input and output
      real(8), intent(in) ::  x_def(:, :), x_ref(:, :), u_c(3), t(3)
      real(8), intent(in) ::  r(3, 3), psi(:)
      real(8)             ::  w
      ! Internal
      integer             :: i, non_zero_nodes, sz(2)
      real(8)             :: zero_tol
      parameter              (zero_tol=1.d-6)
      ! Function start
      sz = shape(x_def)
      w = 0.d0
      non_zero_nodes = 0
      do i = 1, sz(2)
        if (abs(psi(i)) > zero_tol) then
          non_zero_nodes = non_zero_nodes + 1
          w=w+dot_product(t,x_def(:,i)-matmul(r,x_ref(:, i)+u_c))/psi(i)
        end if
      end do
      w = w / non_zero_nodes
      end function
C ******************************************************************** C
C     Compute the linearized warping vector
C ******************************************************************** C
      ! Returns the linearized warping vector.
      ! @ input c: Location of nodes in initial config. relative to 
      !            centroid, size 3xN
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input psi: Warping function for all the nodes, size N
      ! @ input r0: Rotation from referece to initial config., 3x3
      ! @ input r: Rotation from initial to deformed config., 3x3
      ! @ input n_area: Tributary area for each node, size N
      ! @ returns: Linearization of warping amplitude w.r.t x, size 3Nx1
      !
      ! Notes:
      !   - we containts the area weights for each node, then the shear
      !   weights for each node, only use the first set
      pure function calc_lin_w(c, n_sh, psi, r0, r, n_area)
     1              result(w_linear)
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   ::  c(3, n_sh), psi(n_sh), r0(3,3), r(3, 3)
      real(8), intent(in)   ::  n_area(n_sh)
      real(8)               ::  w_linear(3*n_sh)
      ! Internal
      integer               ::  i
      real(8)               ::  bimom, w_tensor(3, 3), c2(3, n_sh)
      real(8)               ::  zero_tol
      parameter                 (zero_tol=1.d-6)
      ! Rotate from the initial config to the reference config
      c2 = matmul(transpose(r0), c)
      ! Assemble the linearized warping vector
      w_linear = 0.d0
      bimom = 0.d0
      do i = 1, n_sh
        if (abs(psi(i)) > zero_tol) then
          ! todo: bimom can be calculated in the function with psi
          bimom = bimom + abs(c2(1, i) * psi(i) * n_area(i))
          w_tensor = 0.d0
          ! Warping only causes axial stress (in 33 direction)
          w_tensor(3, 3) = psi(i) * n_area(i)
          w_linear(3*i-2:3*i) = weight_rotater(w_tensor, r0, r)
        end if
      end do
      ! Normalize by the total bimoment to satisfy equilibrium
      ! B = (M_f,t + M_f,b) * (d_cl / 2)
      d_cl = maxval(c2(2, :))
      bimom = bimom * d_cl
      w_linear = w_linear / bimom      
      end function
C ******************************************************************** C
C     Rotate weight tensors into the deformed configuration
C ******************************************************************** C
      ! Returns the weight vector in the deformed configuration given
      !   a tensor in the reference configuration.
      ! @ input w_ref: Weight tensor in reference config. (3, 3)
      ! @ input r0: Rotation from reference to initial config. (3, 3)
      ! @ input r: Rotation from initial to deformed config. (3, 3)
      ! @ returns: Weight vector normal to cross-section (3, )
      pure function weight_rotater(w_ref, r0, r) result(w_def)
      ! Input and output
      real(8), intent(in) ::  w_ref(3, 3), r0(3, 3), r(3, 3)
      real(8)             ::  w_def(3)
      ! Internal
      real(8)             ::  t0(3), rtot(3, 3)
      ! Function start
      t0 = r0(:, 3)
      rtot = matmul(r, r0)
      w_def = matmul(matmul(rtot, matmul(w_ref, transpose(rtot))), t0)
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
      !   - tf and tw are defined as tf = t_f; tw = E_w / E_f * t_w 
      !     to account for differences in elastic moduli
      !   - The weights are only valid if the interface is elastic
      !   - The flange and web mesh distances are assumed to be constant
      !     between nodes
      pure function calc_weights(p, n_sh, tf_tw_code) result(w)
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   ::  p(3, n_sh)
      integer, intent(in)   ::  tf_tw_code
      real(8)               ::  w(4, n_sh)
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
      integer               ::  i, d_digits, deci, n_digits, c, tf_tw
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
      do i = 1, n_sh
        pt = p(:, i)
        if (abs(pt(1)) < zero_tol) then
          ! Point is on the web
          classification(i) = web_tag
        end if
        if (abs(abs(pt(2)) - y_max) < zero_tol) then
          ! Point is on the flange
          if (classification(i) == web_tag) then
            ! Point is on the web and flange
            classification(i) = joint_tag;
          else if (abs(abs(pt(1)) - x_max) < zero_tol) then
            ! Point is on the corner of the flange
            classification(i) = corner_tag
          else
            ! Just on the flange
            classification(i) = flange_tag
          end if
        end if
      end do
      
      ! Calculate the mesh distance for flange and web nodes
      ! Start search for 2nd greatest values
      x_max_2 = -x_max
      y_max_2 = -y_max
      do i = 1, n_sh
        pt = p(:, i)
        c = classification(i)
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
      ! The thickness values are what is left after we subtract the deci 
      ! value from the front
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
     1                          delta_f, delta_w, classification)
      ! Assign the area weights
      do i = 1, n_sh
        c = classification(i)
        if (c == web_tag) then
          w(1, i) = web_weight
        else if (c == flange_tag) then
          w(1, i) = flange_weight
        else if (c == corner_tag) then
          w(1, i) = corner_weight
        else if (c == joint_tag) then
          w(1, i) = joint_weight
        end if
        ! Add the shear weights
        w(2:4, i) = v_weights(i, :)
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
        v = sgn * (bf - 2. * abs(x1)) * d_cl / 4.
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
        v(1) = (2.*bf*(d+d1)*tf + tw*(d1**2-4.*y1**2)) / (8.*tw)
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
      ! @ returns: The normalized shear resultants, size Nx3
      !   - Column 1: strong axis
      !   - Column 2: weak axis
      !   - Column 3: complementary force for strong axis
      pure function shear_factors(p, n_sh, d, bf, tf, tw, 
     1                            delta_f, delta_w, c_tags) result(v)
      ! Input and output
      integer, intent(in) ::  n_sh, c_tags(n_sh)
      real(8), intent(in) ::  p(3, n_sh), d, bf, tf, tw, 
     1                        delta_f, delta_w
      real(8)             ::  v(n_sh, 3)
      ! Internal variables
      real(8)             ::  d1, tau(2), total_strong, total_weak,
     1                        p_corner(3), tau_2, d_cl
      integer             ::  i, c
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
      do i = 1, n_sh
        c = c_tags(i)
        if (c == web_tag) then
          tau = tw * delta_w * shear_web(p(:, i), d, d1, bf, tf, tw)
          ! Don't consider the complementary shear for the web
          tau_2 = 0.
        else if (c == flange_tag) then
          tau = tf * delta_f * shear_flange(p(:, i), d, d1, bf, tf, tw)
          tau_2 = tf * delta_f * comp_shear_flange(p(:, i), d_cl, bf)
        else if (c == corner_tag) then
          ! Sample at one-half element distance towards the web
          p_corner(1) = -sign(1., p(1, i)) * delta_f / 2.
          p_corner(2:3) = 0.
          p_corner = p(:, i) + p_corner
          tau = tf*delta_f/2.* shear_flange(p_corner, d, d1, bf, tf, tw)
          tau_2 = tf * delta_f / 2.* comp_shear_flange(p_corner,d_cl,bf)
        else
          tau = tf * delta_f * shear_flange(p(:, i), d, d1, bf, tf, tw)
          tau = tau+ tw*delta_w/2.*shear_web(p(:, i), d, d1, bf, tf, tw)
          ! Complementary shear cancels from both sides at the joint
          tau_2 = 0.
        end if
        total_strong = total_strong + tau(1)
        total_weak = total_weak + tau(2)
        v(i, 1:2) = tau
        v(i, 3) = tau_2
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
      real(8), intent(in) ::  w_all(4, n_sh), a_total, r0(3, 3), r(3, 3)
      integer, intent(in) ::  n_sh
      real(8)             ::  u_lin(3, 3 * n_sh)
      ! Internal variables
      real(8)             ::  we(3, 3), rr0(3, 3)
      integer             ::  i   
      
      ! Function start
      do i = 1, n_sh
        we = 0.d0
        ! Weak axis force due to weak axis displacement
        we(1, 1) = w_all(3, i)
        ! Weak axis force in flange due to strong axis displacement
        we(2, 1) = w_all(4, i)
        ! Strong axis force due to strong axis displacement
        we(2, 2) = w_all(2, i)
        ! Axial force due to axial displacement
        we(3, 3) = w_all(1, i) / a_total
        ! Rotate the weights from the reference to the deformed config
        ! w_rotated = r0.t * r.t * we * r * r0
        rr0 = matmul(r, r0)
        u_lin(:, 3*i-2:3*i) = matmul(rr0, matmul(we, transpose(rr0)))
      end do
      end function
C ******************************************************************** C
      end  ! END SUBROUTINE
