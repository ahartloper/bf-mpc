! Modules for the MPC subroutine in the best-fit MPC for beam-to-shell coupling

      include "jacobi_eigenvalue.f90"

      module mpc_modules

      contains
! ******************************************************************** !
  
      pure function skew_sym(v) result(vss)
        ! Returns the 3x3 skew symmetric matrix from axial vector v
        ! @ input v: Axial vector of length 3
        ! @ returns vss: 3x3 skew symmetric matrix
        implicit none
        real(8), intent(in) :: v(3)
        real(8)             :: vss(3, 3)
        ! Specify the columns of the skew_symmetric matrix
        vss(:, 1) = [ 0.d0, v(3), -v(2) ]
        vss(:, 2) = [ -v(3), 0.d0, v(1) ]
        vss(:, 3) = [ v(2), -v(1), 0.d0 ]
      end function

! ******************************************************************** !
  
      pure function left_quat(v0, v) result(ql)
      ! Returns the 4x4 left-side matrix representing {v0, v}
      ! @ input v0: Real part of rotation quaternion
      ! @ input v: Imaginary part of rotation quaternion, length 3
      ! @ return ql: Matrix representing left-hand quaternion 
      !              multiplication
        implicit none
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

! ******************************************************************** !
  
      pure function right_quat(v0, v) result(qr)
      ! Returns the 4x4 right-side matrix representing {v0, v}
      ! @ input v0: Real part of rotation quaternion
      ! @ input v: Imaginary part of rotation quaternion, length 3
      ! @ return qr: Matrix representing right-hand quaternion 
      !              multiplication
        implicit none
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

! ******************************************************************** !
  
      pure function rot_mat(q) result(r)
      ! Returns the 3x3 matrix representing rotation quaternion q
      ! @ input q: Length 4 rotation quaternion
      ! @ returns r: 3x3 rotation matrix
      implicit none
        real(8), intent(in) :: q(4)
        real(8)             :: r(3, 3), one, two
        parameter(one=1.d0, two=2.d0)
        ! Specify the columns
        r(:, 1) = [ one - two * (q(3) ** 2 + q(4) ** 2), &
                    two * (q(2) * q(3) + q(1) * q(4)), &
                    two * (q(2) * q(4) - q(1) * q(3)) ]
        r(:, 2) = [ two * (q(2) * q(3) - q(1) * q(4)), &
                    one - two * (q(2) ** 2 + q(4) ** 2), &
                    two * (q(3) * q(4) + q(1) * q(2)) ]
        r(:, 3) = [ two * (q(2) * q(4) + q(1) * q(3)), &
                    two * (q(3) * q(4) - q(1) * q(2)),&
                    one - two * (q(2) ** 2 + q(3) ** 2) ]
      end function

! ******************************************************************** !
  
      pure function skew_mat(r) result(rr)
      ! Returns the 9x3 skew symmetric matrix from the columns of R
      ! @ input r: 3x3 matrix
      ! @ outputs rr: Vertical stack of the skew symmetric matrices 
      !               formed from the 3 columns of R
      implicit none
        real(8), intent(in) :: r(3, 3)
        real(8)             :: rr(9, 3)
        !
        rr(1:3, :) = skew_sym(r(:, 1))
        rr(4:6, :) = skew_sym(r(:, 2))
        rr(7:9, :) = skew_sym(r(:, 3))
      end function

! ******************************************************************** !
  
      pure function vec2mat_9(v) result(r)
      ! Returns the 3x3 matrix from the vector v
      ! @ input v: Vector of lenght 9
      ! @ returns r: 3x3 matrix
      implicit none
        real(8), intent(in) :: v(9)
        real(8)             :: r(3, 3)
        !
        r(:, 1) = [v(1), v(4), v(7)]
        r(:, 2) = [v(2), v(5), v(8)]
        r(:, 3) = [v(3), v(6), v(9)]
      end function

! ******************************************************************** !
  
      function calc_opquat(b) result(lambda_and_q)
        ! Returns the optimal rotation and associated eigenvalue
        ! @ input b: \mathcal{B} matrix from [2]
        ! @ returns lambda_and_q: First entry is the minimum eigenvalue,
        !                         following entries are the associated 
        !                         eigenvector
        implicit none
        real(8), intent(in) ::  b(4, 4)
        real(8)             ::  be(4, 4), v(4, 4), d(4)
        real(8)             ::  lambda_and_q(5)
        integer             ::  n, it_max, it_num, rot_num
        parameter               (n = 4, it_max = 50)
        !
        be = b  
        call jacobi_eigenvalue(n, b, it_max, v, d, it_num, rot_num)
        lambda_and_q(1) = d(4)
        lambda_and_q(2:5) = v(:, 4)
      end function

! ******************************************************************** !

  ! This function is not required if using the nodal force method
  ! The solution of a 4x4 system is required within (previously done
  ! using LAPACK).

  ! function calc_g(q, n_sh, b, c, lam, r) result(g)
  !   ! Returns the instantaneous rotation matrix
  !   ! @ input q: Optimal rotation quaternion
  !   ! @ input n_sh: Number of shell nodes (N)
  !   ! @ input b: \matcal{B} matrix from [2]
  !   ! @ input c: C matrix from [2]
  !   ! @ input lam: Minimum eigenvalue of B matrix
  !   ! @ input r: Rotation matrix of optimal rotation quaternion
  !   ! @ returns g: Instantaneous rotation matrix
  !   ! Input and output
  !   implicit none
  !   integer, intent(in)   ::  n_sh
  !   real(8), intent(in)   ::  q(4), b(4, 4), c(3, n_sh), r(3, 3)
  !   real(8), intent(in)   ::  lam
  !   real(8)               ::  g(3, 3*n_sh)
  !   ! Internal
  !   integer               ::  i
  !   real(8)               ::  qrr(4, 4), qrr_3(4, 3), xinv(3, 3), id4(4, 4)
  !   real(8)               ::  a_temp(3, n_sh)
  !   ! LAPACK
  !   integer               :: n_ls, nrhs
  !   parameter                (n_ls = 3)
  !   integer               :: lda_ls, ldb_ls
  !   parameter                (lda_ls = 3, ldb_ls = 3)
  !   integer               :: info
  !   integer               :: ipiv(n_ls)
  !   external dgesv
  !   ! Function start
  !   nrhs = 3 * n_sh
  !   ! Calculate I_4x4
  !   id4(:, :) = 0.d0
  !   forall(i = 1:4) id4(i, i) = 1.d0
  !   ! Calculate instantaneous rotation
  !   qrr = right_quat(q(1), q(2:4))
  !   qrr_3 = qrr(:, 2:4)
  !   xinv = matmul(matmul(transpose(qrr_3), (b - lam * id4)), qrr_3)
  !   a_temp = matmul(r, c)
  !   do i = 1, n_sh
  !     g(:, 3*i-2:3*i) = transpose(-skew_sym(a_temp(:, i)))
  !   end do
  !   call dgesv( n_ls, nrhs, xinv, lda_ls, ipiv, g, ldb_ls, info )
  !   g = 4.d0 * g
    
  ! end function

! ******************************************************************** !
  
      pure function calc_q(q, n_sh, g) result(qq)
      implicit none
      ! Returns the linearization of the optimal rotation vector
      ! @ input q: Optimal rotation quaternion, size 4
      ! @ input n_sh: Number of shell nodes (N)
      ! @ input g: Instantaneous rotation matrix, 3x3N
      ! @ returns qq: Linearized rotation matrix
      ! Input and output
      integer, intent(in)   ::  n_sh
      real(8), intent(in)   ::  q(4), g(3, 3*n_sh)
      real(8)               ::  qq(3, 3*n_sh)
      ! Internal
      real(8)               ::  qrr(4, 4), qrr_3(4, 3)
      integer               ::  quatdim
      parameter(quatdim=4)
      ! Function start
      ! Compute the linearized rotation matrix
      qrr = right_quat(q(1), q(2:quatdim))
      qrr_3 = qrr(:, 2:quatdim)
      qq = matmul(qrr_3(2:quatdim, :), g)
      
      end function

! ******************************************************************** !
      
      pure function order_ini(o) result(o2)
      ! Orders o to have the orientation x, y, z as the columns of o2
      ! @ input o: Un-ordered reference configuration
      ! @ returns o2: Ordered reference configuration such that the 
      !               columns of o2 are a right-handed coordinate system
      implicit none
      real(8), intent(in) :: o(3, 3)
      real(8)             :: o2(3, 3), xo(3), yo(3), zo(3), cross_vec(3)
      real(8)             :: test
      ! z is the min eigenvalue, y is the greatest eigenvalue, 
      ! x is the remaining eigenvalue
      zo = o(:, 3)
      yo = o(:, 1)
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

! ******************************************************************** !

      function ini_config(init_pts, n_sh) result(o)
        ! Returns the right-hand ordered reference configuration
        ! @ input init_pts: initial [x, y, z] coordinates of each point
        ! @ input n_sh: Number of shell nodes (N)
        ! @ returns o: Orientation of the reference configuration that is 
        !              a valid right-handed coordinate system.
      ! Input and output
      implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  init_pts(3, n_sh)
        real(8)             ::  o(3, 3)
        ! Internal
        real(8)             ::  ae(3, 3), v(3, 3), d(3)
        integer             ::  ii, n, it_max, it_num, rot_num
        parameter               (n = 3, it_max = 50)
        !
        ae(:, :) = 0.d0
        ! The spread is used to compute the outer product of two vectors
        do ii = 1, n_sh
          ae = ae + spread(init_pts(:, ii), dim=2, ncopies=3) 
     1     * spread(init_pts(:, ii), dim=1, ncopies=3)
        end do
        call jacobi_eigenvalue(n, ae, it_max, v, d, it_num, rot_num)
        ! The eigenvalues are returned in descending order
        o = order_ini(v)
      end function

! ******************************************************************** !
  
      pure function warp_fun(c, n_sh) result(psi)
      ! Returns the warping function evaluated at each node.
      ! @ input c: Location of nodes in reference config., size 3xN
      ! @ input n_sh: Number of shell nodes (N)
      ! @ returns psi: Warping function at each node.
        ! Input and output
        implicit none
        integer, intent(in)   ::  n_sh
        real(8), intent(in)   ::  c(3, n_sh)
        real(8)               ::  psi(n_sh)
        integer               ::  ii
        ! Calculate the warping function, \psi = X*Y
        do ii = 1, n_sh
            psi(ii) = c(1, ii) * c(2, ii)
        end do
      end function

! ******************************************************************** !
  
      pure function warp_amp(x_def, t, psi) result(w)
      ! Returns the value of the averaged warping amplitude.
      ! @ input x_def: pos. of each node in the deformed configuration
      ! @ input t: Orientation of the normal to the cross-section
      ! @ input psi: Warping function at each node
      ! @ returns w: Warping amplitude
        ! Input and output
        implicit none
        real(8), intent(in) ::  x_def(:, :), t(3)
        real(8), intent(in) ::  psi(:)
        real(8)             ::  w
        ! Internal
        integer             :: ii, sz(2)
        real(8)             :: psi2
        ! Function start
        sz = shape(x_def)
        w = 0.d0
        psi2 = 0.d0
        do ii = 1, sz(2)
          w = w + psi(ii) * t(1)*x_def(1, ii) + t(2)*x_def(2, ii) 
     1       + t(3)*x_def(3, ii)
          psi2 = psi2 + psi(ii)*psi(ii)
        end do
        w = w / psi2
      end function

! ******************************************************************** !
  
      pure function extract_rotation(q) result(phi)
      ! Returns the rotation vector extracted from the quaternion
      ! @ input q: Rotation quaternion
      ! @ returns phi: Euler angle representation of q
      !
      ! From Abaqus Theory Guide 1.3.1
        ! Input and output
        implicit none
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

! ******************************************************************** !
  
      pure function shear_flange(xi, d, d1, bf, tf, tw) result(v)
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
        ! Input and output
        implicit none
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

  ! ******************************************************************** !
  
      pure function comp_shear_flange(xi, d_cl, bf) result(v)
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
        ! Input and output
        implicit none
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

! ******************************************************************** !
  
      pure function shear_web(xi, d, d1, bf, tf, tw) result(v)
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
        ! Input and output
        implicit none
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

! ******************************************************************** !
  
      pure function shear_factors(p, n_sh, sec_props, mesh_sizes, 
     1                            c_tags) result(v)
        ! Returns the shear factor resultant for one flange and web
        ! @ input p: Shell nodes in the reference configuration, 3xN
        ! @ input n_sh: Number of shell nodes (N)
        ! @ input d: Section depth
        ! @ input bf: Section width
        ! @ input tf: Flange thickness
        ! @ input tw: Web thickness
        ! @ input c_tags: Node classification tag, size N
        ! @ input etang: Tangent modulus, size N
        ! @ returns: The normalized shear resultants, size Nx4
        !   - Column 1: strong axis
        !   - Column 2: weak axis
        !   - Column 3: complementary force for strong axis
        ! Input and output
        implicit none
        integer, intent(in) ::  n_sh, c_tags(n_sh)
        real(8), intent(in) ::  p(3, n_sh), sec_props(4), mesh_sizes(2)
        real(8)             ::  v(n_sh, 3)
        ! Internal variables
        real(8)             ::  d1, tau(2), total_strong, total_weak
        real(8)             ::  p_corner(3), tau_2, d_cl, tau_3, 
       1                        tau_temp(2), p_temp(3)
        real(8)             ::  delta_f, delta_w, d, bf, tf, tw
        integer             ::  ii, c
        ! Classification tags
        integer             ::  web_tag, flange_tag, corner_tag, 
       1                        joint_tag
        parameter              (web_tag=1, flange_tag=2, corner_tag=3,
       1                        joint_tag=4)
        ! Function start
        d = sec_props(1)
        bf = sec_props(2)
        tf = sec_props(3)
        tw = sec_props(4)
        delta_f = mesh_sizes(1)
        delta_w = mesh_sizes(2)
        ! The first value of tau is for the strong axis, the second is for
        ! the weak axis
        ! tau_2 is the complementary shear in the flange
        ! tau_3 is the normal direction of the shear stress resultant
        d1 = d - 2. * tf
        d_cl = d - tf
        total_strong = 0.
        total_weak = 0.
        do ii = 1, n_sh
          c = c_tags(ii)
          if (c == web_tag) then
            tau = tw * delta_w * shear_web(p(:, ii),d,d1,bf,tf,tw)
            ! No complementary shear for the web
            tau_2 = 0.
            
          else if (c == flange_tag) then
            tau = tf * delta_f * shear_flange(p(:, ii),d,d1,bf,tf,tw)
            tau_2 = tf * delta_f * comp_shear_flange(p(:, ii),d_cl,bf)
            
          else if (c == corner_tag) then
            ! Sample at one-half element distance towards the web
            p_corner(1) = -sign(1.d0, p(1, ii)) * delta_f / 2.
            p_corner(2:3) = 0.
            p_corner = p(:, ii) + p_corner
            tau = tf*delta_f/2.* shear_flange(p_corner,d,d1,bf,tf,tw)
            tau_2 = tf * delta_f / 2.* comp_shear_flange(p_corner,d_cl,bf)
            
          else  ! joint node
            tau=tf*delta_f*shear_flange(p(:, ii),d,d1,bf,tf,tw)
            
            tau = tau+tw*delta_w/2.*shear_web(p(:, ii),d,d1,bf,tf,tw)
            ! Complementary shear cancels from both sides at the joint
            tau_2 = 0.
          end if
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
  
! ******************************************************************** !

      pure function calc_disp_lin(n_sh, w_all, a_total, r0, r) result(u_lin)
      ! Returns the matrix relating centroidal disp. to shell disps.
      ! @ input w_all: Area and shear weights, size 4xN
      ! @ input a_total: Total area
      ! @ input r0: Initial orientation, size 3x3
      ! @ input r: Rotation from initial to deformed, size 3x3
      ! @ returns: The displacement factor for each direction for each 
      !            node, size 3xN
        ! Input and output
        implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  w_all(4, n_sh), a_total, r0(3, 3),
       1                        r(3, 3)
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
          ! Weak axis nodal force due to weak axis shear force
          w_bar = 0.d0
          w_bar(1, 3) = w_all(3, ii)
          w_bar(3, 1) = w_all(3, ii)
          s_bar(:, 1) = matmul(w_bar, t_bar)
          ! Weak axis nodal force in flange due to strong axis shear force
          w_bar = 0.d0
          w_bar(3, 1) = w_all(4, ii)
          w_bar(1, 3) = w_all(4, ii)
          ! Strong axis nodal force due to strong axis shear force
          w_bar(3, 2) = w_all(2, ii)
          w_bar(2, 3) = w_all(2, ii)
          s_bar(:, 2) = matmul(w_bar, t_bar)
          ! Axialnodal force due to axial force
          w_bar = 0.d0
          w_bar(3, 3) = w_all(1, ii) / a_total
          s_bar(:, 3) = matmul(w_bar, t_bar)
          ! Rotate the nodal forces from the reference to the deformed config
          s = matmul(rr0, matmul(s_bar, rr0_T))
          u_lin(1:3, 3*ii-2:3*ii) = transpose(s)
        end do
      end function

! ******************************************************************** !

      function torquer(n_sh, p_ref, areas, d_cl, bf) result(tq)
        ! Returns torque on flange only
        ! @input n_sh: Number of shell nodes on the interface
        ! @input xyz: Coords of shell nodes in reference config.
        ! @input areas: Node areas.
        ! @input d_cl: Section centerline depth.
        ! @input bf: Section flange width.
        ! @returns: Nodal forces due to a unit torque on the cross-section.
        ! Input and output
        implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  p_ref(3, n_sh), areas(n_sh), d_cl, bf
        real(8)             ::  tq(2, n_sh)
        real(8)             ::  tor, ff
        integer             ::  ii
        real(8)             ::  tol
        parameter(tol=1.d-3)
        ! Function    
        tor = 0.
        tq = 0.d0
        do ii = 1, n_sh
          if ( abs(abs(p_ref(2, ii)) - (d_cl / 2.d0)) < tol ) then
            tq(1:2, ii) = areas(ii) * [ -p_ref(2, ii), p_ref(1, ii) ]
            tor = tor + areas(ii) * norm2(p_ref(1:2, ii)) ** 2
          end if
        end do
        tq = tq / tor
      end function

! ******************************************************************** !

      pure function lin_rot(n_sh, r0, r, rw) result(linr)
        ! Returns the linearization of the rotation constraints.
        ! @input n_sh: Number of shell nodes on the interface
        ! @input r0: Rotation matrix from reference to initial configs
        ! @input r: Rotation matrix from initial to current configs
        ! @input rw: Weights for the moments / torque
        ! @returns: Matrix for the linearization of the rotations
        !
        ! Notes:
        !   - rw ordering:
        !     - rw(1, :) = axial due to moment about x
        !     - rw(2, :) = axial due to moment about y
        !     - rw(3, :) = shear in x due to moment about z (torsion)
        !     - rw(4, :) = shear in y due to moment about z (torsion)
        implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  r(3, 3), r0(3, 3), rw(4, n_sh)
        real(8)             ::  linr(3, 3*n_sh)
        integer             ::  ii
        real(8)             ::  w_bar(3, 3), s_bar(3, 3), s(3, 3), 
       1                        t_bar(3), rr0(3, 3), rr0_T(3, 3)
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
      end function
  
! ******************************************************************** !

      pure function normalized_coords(n_sh, xy, d, bf) result(xy_bar)
        ! Returns the normalized x,y-coordinates of all the nodes.
        ! @input n_sh: Number of shell elements.
        ! @input xy_bar: Node reference x,y-cooridinates
        ! @input d: Section depth.
        ! @input bf: Web width.
        ! @returns: Normlized node reference x,y-cooridinates.
        implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  xy(2, n_sh), d, bf
        real(8)             ::  xy_bar(2, n_sh)
        integer             ::  i
        ! Function start
        do i = 1, n_sh
          xy_bar(1, i) = xy(1, i) / (bf / 2.)
          xy_bar(2, i) = xy(2, i) / (d / 2.)
        end do
      end function

! ******************************************************************** !

      pure function compute_section_props(n_sh, xy_bar, da) result(sp)
        ! Returns the area and normalized moments of inertia that are needed for the section stiffness.
        ! @input n_sh: Number of shell elements.
        ! @input xy_bar: Normlized node reference x,y-cooridinates.
        ! @input da: Node tributary areas.
        implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  xy_bar(2, n_sh), da(n_sh)
        real(8)             ::  sp(4)
        integer             ::  i
        ! Function start
        sp = 0.d0
        do i = 1, n_sh
          sp(1) = sp(1) + da(i)
          sp(2) = sp(2) + xy_bar(2, i)**2 * da(i)
          sp(3) = sp(3) + xy_bar(1, i)**2 * da(i)
          sp(4) = sp(4) + xy_bar(1, i)**2 * xy_bar(2, i)**2 * da(i)
        end do
      end function

! ******************************************************************** !

      pure function section_tangent(n_sh, xy_bar, da, etang, section_props) result(ksec)
        ! Returns the section tangent stiffness matrix.
        ! @input n_sh: Number of shell elements.
        ! @input xy_bar: Normlized node reference x,y-cooridinates.
        ! @input da: Node tributary areas.
        ! @input etang: Tangent modulus for each node.
        ! @input section_props: Normalized [A,Ix,Iy,Iw].
        ! @returns: Assembled tangent stiffness matrix.
        !
        ! Notes:
        !   - Based on (6.22) on pg. 277 from Chen and Atsuta 1977.
        implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  xy_bar(2, n_sh),da(n_sh),etang(n_sh),
       1                        section_props(4)
        real(8)             ::  ksec(4, 4), ks(4, 4)
        real(8)             ::  x, y, a_bar, ix_bar, iy_bar, iw_bar
        integer             ::  i
        ! Function start
        a_bar = section_props(1)
        ix_bar = section_props(2)
        iy_bar = section_props(3)
        iw_bar = section_props(4)
        ksec = 0.d0
        do i = 1, n_sh
          x = xy_bar(1, i)
          y = xy_bar(2, i)
          ! Column 1
          ks(1, 1) = etang(i) * da(i) / a_bar
          ks(2, 1) = etang(i) * y * da(i) / ix_bar
          ks(3, 1) = -etang(i) * x * da(i) / iy_bar
          ks(4, 1) = etang(i) * x * y * da(i) / iw_bar
          ! Column 2
          ks(1, 2) = etang(i) * y * da(i) / a_bar
          ks(2, 2) = etang(i) * y**2 * da(i) / ix_bar
          ks(3, 2) = -etang(i) * x * y * da(i) / iy_bar
          ks(4, 2) = etang(i) * x * y**2 * da(i) / iw_bar
          ! Column 3
          ks(1, 3) = -etang(i) * x * da(i) / a_bar
          ks(2, 3) = -etang(i) * x * y * da(i) / ix_bar
          ks(3, 3) = etang(i) * x**2 * da(i) / iy_bar
          ks(4, 3) = -etang(i) * x**2 * y * da(i) / iw_bar
          ! Column 4
          ks(1, 4) = etang(i) * x * y * da(i) / a_bar
          ks(2, 4) = etang(i) * x * y**2 * da(i) / ix_bar
          ks(3, 4) = -etang(i) * x**2 * y * da(i) / iy_bar
          ks(4, 4) = etang(i) * x**2 * y**2 * da(i) / iw_bar
          ksec = ksec + ks
        end do
      end function
  
! ******************************************************************** !

      pure function calc_section_deform(ksec, sec_props, d, bf) 
     1              result(usec)
        ! Returns the generalized section deformations for unit loads in each DOF.
        ! @input ksec: Normalized section stiffness matrix.
        ! @input sec_props: "Regular" section properties
        ! @input d: Section depth
        ! @input bf: Section width
        ! @returns: Section deformations, column i contains the defomations due to a force in DOF i.
        ! 
        ! Notes:
        !   - Deformations: [centroid axial strain, strong axis curvature, weak axis curvature, torsion warping curvature]
        !   - sec_props: [A, Ix, Iy, Iw]
        implicit none
        real(8), intent(in) ::  ksec(4, 4), sec_props(4), d, bf
        real(8)             ::  usec(4, 4)
        real(8)             ::  kinv(4, 4), f(4), 
     1                          normalization_factor(4)
        integer             ::  i
        ! Function start
        normalization_factor(1) = sec_props(1)
        normalization_factor(2) = sec_props(2) / (d / 2.d0)
        normalization_factor(3) = sec_props(3) / (bf / 2.d0)
        normalization_factor(4) = sec_props(4) / (d * bf / 4.d0)
        ! Here we assume that the section is elastic -> section stiffness is identity matrix
        ! Therefore, u_i = f_i
        do i = 1, 4
          f = 0.d0
          f(i) = 1.d0 /  normalization_factor(i)
          usec(1:4, i) = f
        end do
      end function
  
! ******************************************************************** !

      pure function nodal_forces(n_sh, xy_bar, da, etang, usec) result(f_nodes)
        ! Returns the axial force at each node for the deformation in each DOF.
        ! @input n_sh: Number of shell elements.
        ! @input xy_bar: Normlized node reference x,y-cooridinates.
        ! @input da: Node tributary areas.
        ! @input etang: Tangent modulus for each node.
        ! @input usec: Section deformations, column i contains the defomations due to a force in DOF i.
        ! @returns: Nodal forces, column j is for node j, row i is for DOF i.
        implicit none
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  xy_bar(2, n_sh), da(n_sh), etang(n_sh), 
     1                          usec(4, 4)
        real(8)             ::  f_nodes(4, n_sh)
        integer             ::  i, j
        real(8)             ::  eda, x, y
        ! Function start
        do j = 1, n_sh
          x = xy_bar(1, j)
          y = xy_bar(2, j)
          eda = da(j) * etang(j)
          do i = 1, 4
            f_nodes(i, j) = eda * (usec(1, i) + y*usec(2, i) 
     1                      - x*usec(3, i) + x*y*usec(4, i))
          end do
        end do
      end function
  
! ******************************************************************** !

      pure function compute_node_areas(n_sh, xyz, sec_props, 
     1                                classification, mesh_sizes) 
     2                                result(areas)
        ! Returns the area for each node.
        ! @input n_sh: Number of shell elements.
        ! @input xyz: Node coordinates in reference config.
        ! @input sec_props: Cross-section dimensions.
        ! @input classification: Node classification tags.
        ! @input: Mesh sizes for the flange and web, [delta_f, delta_w].
        ! @returns: Tributary area for each node.
        implicit none
        integer, intent(in) ::  n_sh, classification(n_sh)
        real(8), intent(in) ::  xyz(3, n_sh), sec_props(4), mesh_sizes(2)
        real(8)             ::  areas(n_sh)
        integer             ::  i
        real(8)             ::  web_area, flange_area, tol, d_cl, bf
        integer             ::  web_tag, flange_tag, corner_tag, 
     1                          joint_tag
        parameter               (web_tag=1, flange_tag=2, corner_tag=3, 
     1                           joint_tag=4)
        parameter(tol=1.d-3)
        ! Function start
        d_cl = sec_props(1) - sec_props(3)
        bf = sec_props(2)
        flange_area = sec_props(3) * mesh_sizes(1)
        web_area = sec_props(4) * mesh_sizes(2)
        do i = 1, n_sh
          if (classification(i) == web_tag) then
            areas(i) = web_area
          else if(classification(i) == flange_tag) then
            areas(i) = flange_area
          else if(classification(i) == joint_tag) then
            areas(i) = flange_area + 0.5d0 * web_area
          else if(classification(i) == corner_tag) then
            areas(i) = 0.5d0 * flange_area
          end if
        end do
      end function

! ******************************************************************** !

      pure function compute_mesh_size(n_sh, p, sec_props, 
     1                                classification) result(mesh_sizes)
        implicit none
        ! Returns the mesh size for the flange and web on the interface.
        ! @input n_sh: Number of shell elements.
        ! @input p: Node coordinates in reference config.
        ! @input sec_props: Cross-section dimensions.
        ! @input classification: Node classification tags.
        ! @returns: Mesh sizes for the flange and web, [delta_f, delta_w].
        integer, intent(in) ::  n_sh, classification(n_sh)
        real(8), intent(in) ::  p(3, n_sh), sec_props(4)
        real(8)             ::  mesh_sizes(2)
        integer             ::  i, c
        integer             ::  web_tag,flange_tag,corner_tag,joint_tag
        parameter               (web_tag=1, flange_tag=2, corner_tag=3,
     1                           joint_tag=4)
        real(8)             ::  x_max,x_max_2,y_max,y_max_2,tol,pt(3)
        parameter(tol=1.d-3)
        ! Function start
        x_max = sec_props(2) / 2.
        y_max = (sec_props(1) - sec_props(3)) / 2.
        ! get mesh size
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
          if ((c==flange_tag .or. c==joint_tag) .and. pt(2)>0.d0) then
            ! Set the 2nd largest x value
            if (x_max_2 < pt(1) .and. pt(1) < x_max) then
              x_max_2 = pt(1)
            end if
          end if
        end do
        ! Flange then web
        mesh_sizes(1) = x_max - x_max_2
        mesh_sizes(2) = y_max - y_max_2
      end function

! ******************************************************************** !
  
      pure function classify_nodes(n_sh, xyz, d_cl, bf) 
     1              result(classification)
        implicit none
        ! Returns classification for each node.
        ! @input xyz: Nodes in reference configuration.
        ! @input d_cl: Centerline depth.
        ! @input bf: Web width.
        ! @returns: The classification for each node.
        integer, intent(in) ::  n_sh
        real(8), intent(in) ::  xyz(3, n_sh), d_cl, bf
        integer             ::  classification(n_sh)
        integer             ::  i
        integer             ::  web_tag,flange_tag,corner_tag,joint_tag
        parameter               (web_tag=1, flange_tag=2, corner_tag=3,
     1                           joint_tag=4)
        real(8)             ::  tol
        parameter(tol=1.d-3)
        ! Function start
        do i = 1, n_sh
          if (abs(abs(xyz(2, i)) - d_cl / 2.) < tol) then
            ! Flange / joint nodes
            if (abs(xyz(1, i)) < tol) then
              ! Joint node
              classification(i) = joint_tag
            else if (abs(abs(xyz(1, i)) - bf / 2.) < tol) then
              ! End node
              classification(i) = corner_tag
            else
              classification(i) = flange_tag
            end if
          else
            ! Web node
            classification(i) = web_tag
          end if
        end do

      end function

! ******************************************************************** !

end module mpc_modules
