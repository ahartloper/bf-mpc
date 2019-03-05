      SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,
     * LMPC,KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
C
      !INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION UE(MDOF),A(MDOF,MDOF,N),JDOF(MDOF,N),X(6,N),
     * U(MAXDOF,N),UINIT(MAXDOF,N),TIME(2),TEMP(NT,N),
     * FIELD(NF,NT,N),LTRAN(N),TRAN(3,3,N)
C
C *************************************************************************** C
C Beam to shell coupling with warping
C A Hartloper, EPFL
C *************************************************************************** C
      ! Internal variables
      ! Defines the spatial DOF and size of quaternions
      INTEGER               :: NDIM, QUATDIM
      PARAMETER                (NDIM = 3, QUATDIM = 4)
      ! Number of shell nodes
      INTEGER               :: n_shell
      REAL(8)               :: n_reciprocal
      ! Loop counter
      INTEGER               :: i, i_sh
      ! Warping amplitude
      REAL(8)               :: w_amp
      ! Centroids (deformed and reference)
      REAL(8)               :: XC(NDIM), XC0(NDIM), u_cent(NDIM)
      ! First shell node (deformed and translated reference)
      REAL(8)               :: x_s1(NDIM), x_s1_ref(NDIM)
      ! Rotation matrices (to get deformed, reference)
      REAL(8)               :: R_mat(NDIM, NDIM), R_ref(NDIM, NDIM)
      ! Orientation normal to beam plane
      REAL(8)               :: t_def(NDIM)
      ! Arrays that depend on the number of nodes
      REAL(8), ALLOCATABLE  :: C_mat(:, :), D_mat(:, :), G_mat(:, :), 
     1                         Q_mat(:, :), 
     2                         X_shell(:, :), U_shell(:, :), 
     3                         w_lin(:), X_ref(:, :), psi_all(:)
      ! For the rotation quaternion and computations
      REAL(8)               :: B_mat(QUATDIM, QUATDIM), 
     1                         beta(QUATDIM, QUATDIM)
      REAL(8)               :: quat(QUATDIM)
      REAL(8)               :: lambda
      REAL(8)               :: lam_q(QUATDIM+1)
      ! Numerical parameters
      REAL(8)               :: ONE, TWO, ZERO
      PARAMETER                (ONE=1.0D0, TWO=2.0D0, ZERO=0.D0)
      ! Allocate variable size arrays
      n_shell = N - 1
      ALLOCATE(X_shell(NDIM, n_shell))
      ALLOCATE(U_shell(NDIM, n_shell))
      ALLOCATE(X_ref(NDIM, n_shell))
      ALLOCATE(C_mat(NDIM, n_shell))
      ALLOCATE(D_mat(NDIM, n_shell))
      ALLOCATE(G_mat(NDIM, NDIM * n_shell))
      ALLOCATE(Q_mat(NDIM, NDIM * n_shell))
      ALLOCATE(psi_all(n_shell))
      ALLOCATE(w_lin(NDIM * n_shell))
C *************************************************************************** C
C Routine definition
C *************************************************************************** C      
      ! Calculate the centroid locations
      ! The first node is the beam node, the shell nodes follow
      X_shell = X(1:NDIM, 2:N)
      U_shell = U(1:NDIM, 2:N)
      XC0(:) = ZERO
      XC(:) = ZERO
      DO i = 1, n_shell
        XC0(:) = XC0(:) + X_shell(:, i)
        XC(:) = XC(:) + U_shell(:, i)
     1  + X_shell(:, i)
      END DO
      n_reciprocal = ONE / n_shell
      XC0 = XC0 * n_reciprocal
      XC = XC * n_reciprocal
      u_cent = XC - XC0
            
      ! Compute the reference configuration
      DO i = 1, n_shell
        X_ref(:, i) = X_shell(:, i) - XC0
      END DO
      R_ref = refConfig(X_ref)
      
      ! Compute the optimal rotation quaternion
      B_mat(:, :) = ZERO
      DO i = 1, n_shell
        C_mat(:, i) = MATMUL(TRANSPOSE(R_ref), X_shell(:, i) - XC0)
        D_mat(:, i) = X_shell(:, i) + U_shell(:, i) - XC
        beta = leftQuat(ZERO, D_mat(:, i)) - rightQuat(ZERO,C_mat(:, i))
        B_mat = B_mat + MATMUL(TRANSPOSE(beta), beta)
      END DO
      lam_q = getRotQuat(B_mat)
      lambda = lam_q(1)
      quat = lam_q(2:5)
      ! Make the first entry always positive
      quat = SIGN(ONE, quat(1)) * quat  
      
      ! Compute the linearized rotation matrix
      R_mat = rotMat(quat)
      G_mat = getG(quat, B_mat, C_mat, lambda, R_mat)
      Q_mat = getLinR(quat, G_mat)
      
      ! Calculate the warping amplitude
      psi_all = warp_fun(C_mat)
      t_def(:) = MATMUL(R_mat, R_ref(:, 3))
      w_amp = warp_amp(X_shell, X_shell + U_shell, u_cent, t_def, R_mat,
     1                 psi_all)
      ! Calculate the linearized warping
      w_lin = lin_warp(X_shell + U_shell, R_ref(:, 3), psi_all, R_mat, 
     1                 G_mat)
      
      ! Assign the A submatrix for the beam node and set active DOF
      FORALL(i = 1:7) A(i, i, 1) = ONE
      FORALL(i = 1:7) JDOF(i, 1) = i
      
      ! Assign the A submatrices for the shell nodes
      DO i = 2, N
        ! Index corresponding to the shell nodes
        i_sh = i - 1
        ! Displacement constraints
        A(1, 1, i) = -n_reciprocal
        A(2, 2, i) = -n_reciprocal
        A(3, 3, i) = -n_reciprocal
        ! Rotation constraints
        A(4, 1:3, i) = -Q_mat(1, 3*i_sh-2:3*i_sh)
        A(5, 1:3, i) = -Q_mat(2, 3*i_sh-2:3*i_sh)
        A(6, 1:3, i) = -Q_mat(3, 3*i_sh-2:3*i_sh)
        ! Warping constraint
        A(7, 1:3, i) = -w_lin(3*i_sh-2:3*i_sh) 
        ! Set the active DOF in the shell elements
        JDOF(1, i) = 1
        JDOF(2, i) = 2
        JDOF(3, i) = 3
      END DO
            
      ! Update the node DOF exactly
      UE(1:3) = u_cent
      UE(4:6) = extractRotation(quat)
      UE(7) = w_amp
      RETURN
      
      CONTAINS
C *************************************************************************** C
C     Skew symmetric matrix from vector
C *************************************************************************** C
      PURE FUNCTION skewSym(v) result(Vss)
      ! Returns the 3x3 skew symmetric matrix from axial vector v
      ! @ input v: Axial vector of length 3
      ! @ returns Vss: 3x3 skew symmetric matrix
      REAL(8), intent(in) :: v(3)
      REAL(8)             :: Vss(3, 3)
      ! Specify the columns of the skewsymmetric matrix (transpose of the rows)
      Vss(:, 1) = (/ 0.D0, v(3), -v(2) /)
      Vss(:, 2) = (/ -v(3), 0.D0, v(1) /)
      Vss(:, 3) = (/ v(2), -v(1), 0.D0 /)
      END FUNCTION
C *************************************************************************** C
C     Left quaternion representation from 3 element vector
C *************************************************************************** C
      ! Returns the 4x4 left-side matrix representing {v0, v}
      ! @ input v0: Real part of rotation quaternion
      ! @ input v: Imaginary part of rotation quaternion, length 3
      ! @ return Ql: Matrix representing left-hand quaternion multiplication
      PURE FUNCTION leftQuat(v0, v) result(Ql)
      REAL(8), intent(in) :: v0
      REAL(8), intent(in) :: v(3)
      REAL(8)             :: Ql(4, 4)
      INTEGER             :: i
      !
      Ql(:, 1) = (/ v0, v(:) /)
      Ql(1, 2:4) = -v(:)
      Ql(2:4, 2:4) = skewSym(v)
      FORALL(i = 1:3) Ql(1 + i, 1 + i) = Ql(1 + i, 1 + i) + v0
      END FUNCTION
C *************************************************************************** C
C     Right quaternion representation from 3 element vector
C *************************************************************************** C
      ! Returns the 4x4 right-side matrix representing {v0, v}
      ! @ input v0: Real part of rotation quaternion
      ! @ input v: Imaginary part of rotation quaternion, length 3
      ! @ return Qr: Matrix representing right-hand quaternion multiplication
      PURE FUNCTION rightQuat(v0, v) result(Qr)
      REAL(8), intent(in) :: v0
      REAL(8), intent(in) :: v(3)
      REAL(8)             :: Qr(4, 4)
      INTEGER             :: i
      !
      Qr(:, 1) = (/ v0, v(:) /)
      Qr(1, 2:4) = -v(:)
      Qr(2:4, 2:4) = -skewSym(v)   
      FORALL(i = 1:3) Qr(1 + i, 1 + i) = Qr(1 + i, 1 + i) + v0   
      END FUNCTION
C *************************************************************************** C
C     Rotation matrix from quaternion
C *************************************************************************** C
      ! Returns the 3x3 matrix representing rotation quaternion q
      ! @ input q: Length 4 rotation quaternion
      ! @ returns R: 3x3 rotation matrix
      PURE FUNCTION rotMat(q) result(R)
      REAL(8), intent(in) :: q(4)
      REAL(8)             :: R(3, 3)
      ! Specify the columns
      R(:, 1) = (/ ONE - TWO * (q(3) ** 2 + q(4) ** 2), 
     1 TWO * (q(2) * q(3) + q(1) * q(4)), 
     2 TWO * (q(2) * q(4) - q(1) * q(3)) /)
      R(:, 2) = (/ TWO * (q(2) * q(3) - q(1) * q(4)),
     1 ONE - TWO * (q(2) ** 2 + q(4) ** 2),
     2 TWO * (q(3) * q(4) + q(1) * q(2)) /)
      R(:, 3) = (/ TWO * (q(2) * q(4) + q(1) * q(3)),
     1 TWO * (q(3) * q(4) - q(1) * q(2)),
     2 ONE - TWO * (q(2) ** 2 + q(3) ** 2) /)
      END FUNCTION
C *************************************************************************** C
C     Stack of skew symmetric matrix from columns of matrix
C *************************************************************************** C
      PURE FUNCTION skewMat(R) result(RR)
      ! Returns the 9x3 skew symmetric matrix from the columns of R
      ! @ input R: 3x3 matrix
      ! @ outputs RR: Vertical stack of the skew symmetric matrices formed from
      !               the 3 columns of R
      REAL(8), intent(in) :: R(3, 3)
      REAL(8)             :: RR(9, 3)
      !
      RR(1:3, :) = skewSym(R(:, 1))
      RR(4:6, :) = skewSym(R(:, 2))
      RR(7:9, :) = skewSym(R(:, 3))
      END FUNCTION
C *************************************************************************** C
C     Create a 3x3 matrix out of vector of length 9
C *************************************************************************** C
      PURE FUNCTION vec2mat_9(v) result(R)
      ! Returns the 3x3 matrix from the vector v
      ! @ input v: Vector of lenght 9
      ! @ returns R: 3x3 matrix
      REAL(8), intent(in) :: v(9)
      REAL(8)             :: R(3, 3)
      !
      R(:, 1) = v(1:3)
      R(:, 2) = v(4:6)
      R(:, 3) = v(7:9)
      END FUNCTION
C *************************************************************************** C
C     Compute optimal rotation quaternion
C *************************************************************************** C
      FUNCTION getRotQuat(B) result(lambda_and_q)
      ! Returns the optimal rotation and associated eigenvalue
      ! @ input B: \mathcal{B} matrix from Mostafa and Sivaselvan
      ! @ returns lambda_and_q: First entry is the minimum eigenvalue,
      !                         following entries are the associated eigenvector
      REAL(8), intent(in) ::  B(4, 4)
      REAL(8)             ::  BE(4, 4)
      REAL(8)             ::  lambda_and_q(5)
      ! For LAPACK
      INTEGER             ::  NE, NSELECT
      PARAMETER               (NE = 4, NSELECT = 1)
      INTEGER             ::  LDA, LDZ
      PARAMETER               (LDA = NE, LDZ = NE)
      INTEGER             ::  INFO, LWORK, LIWORK, IL, IU, M
      REAL(8)             ::  ABSTOL, VL, VU
      INTEGER :: LWMAX
      PARAMETER(LWMAX = 1000)
      INTEGER             ::  ISUPPZ(NE), IWORK(LWMAX)
      REAL(8)             ::  AE(LDA, NE), W(NE), Z(LDZ, NSELECT), 
     1                        WORK(LWMAX)
      EXTERNAL DSYEVR
      !
      ABSTOL = -1.0
      IL = 1
      IU = NSELECT
      ! LWORK and LIWORK are set based on the query to the optimal 
      ! workspace during testing
      LWORK = 104
      LIWORK = 40
      BE = B  ! We need B for later, so don't overwrite
      CALL DSYEVR('Vectors', 'Indices', 'Upper', NE, BE, LDA, VL, VU,IL,
     1            IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     2            LIWORK, INFO) 
      lambda_and_q(1) = W(1)
      lambda_and_q(2:5) = Z(:, 1)
      END FUNCTION
C *************************************************************************** C
C     Compute G matrix
C *************************************************************************** C
      FUNCTION getG(q, B, C, lam, R) result(G)
      ! Returns the instantaneous rotation matrix from Mostafa and Sivaselvan
      ! @ input q: Optimal rotation quaternion
      ! @ input B: \matcal{B} matrix from Mostafa and Sivaselvan
      ! @ input C: C matrix from Mostafa and Sivaselvan
      ! @ input lam: Minimum eigenvalue of B matrix
      ! @ input R: Rotation matrix of optimal rotation quaternion
      ! @ returns G: Instantaneous rotation matrix
      ! Input and output
      REAL(8), intent(in)   ::  q(4), B(:, :), C(:, :), R(3, 3)
      REAL(8), intent(in)   ::  lam
      REAL(8), ALLOCATABLE  ::  G(:, :)
      ! Internal
      INTEGER               ::  sz(2)
      REAL(8)               ::  qrr(4, 4), qrr_3(4, 3), Xinv(3, 3), 
     1                          ID4(4, 4)
      REAL(8), ALLOCATABLE  ::  A_temp(:, :)
      ! LAPACK
      INTEGER               :: N_LS, NRHS
      PARAMETER                (N_LS = 3)
      INTEGER               :: LDA_LS, LDB_LS
      PARAMETER                (LDA_LS = 3, LDB_LS = 3)
      INTEGER               :: INFO
      INTEGER               :: IPIV(N_LS)
      EXTERNAL DGESV
      ! Assign array sizes
      sz = SHAPE(C)
      ALLOCATE(A_temp(3, sz(2)))
      ALLOCATE(G(3, 3 * sz(2)))
      NRHS = 3 * sz(2)
      ! Function start
      ! Calculate I_4x4
      ID4(:, :) = 0.D0
      FORALL(i = 1:4) ID4(i, i) = 1.D0
      qrr = rightQuat(q(1), q(2:4))
      qrr_3 = qrr(:, 2:4)
      Xinv = MATMUL(MATMUL(TRANSPOSE(qrr_3), (B - lam * ID4)), qrr_3)
      A_temp = MATMUL(R, C)
      DO i = 1, sz(2)
        G(:, 3*i-2:3*i) = TRANSPOSE(-skewSym(A_temp(:, i)))
      END DO
      CALL DGESV( N_LS, NRHS, Xinv, LDA_LS, IPIV, G, LDB_LS, INFO )
      G = 4.D0 * G
      
      END FUNCTION
C *************************************************************************** C
C     Compute linearized rotation matrix
C *************************************************************************** C
      PURE FUNCTION getLinR(q, G) result(QQ)
      ! Returns the linearization of the optimal rotation matrix
      ! @ input q: Optimal rotation quaternion
      ! @ input G: Instantaneous rotation matrix
      ! @ returns QQ: Linearization of rotation matrix R
      ! Input and output
      REAL(8), intent(in)   ::  q(4), G(:, :)
      REAL(8), ALLOCATABLE  ::  QQ(:, :)
      ! Internal
      INTEGER               ::  sz(2)
      REAL(8)               ::  qrr(4, 4), qrr_3(4, 3)      
      ! Assign array sizes
      sz = SHAPE(G)
      ALLOCATE(QQ(3, sz(2)))
      ! Function start
      ! Compute the linearized rotation matrix
      qrr = rightQuat(q(1), q(2:QUATDIM))
      qrr_3 = qrr(:, 2:QUATDIM)
      QQ = MATMUL(qrr_3(2:QUATDIM, :), G)
      
      END FUNCTION
C *************************************************************************** C
C     Order reference configuration
C *************************************************************************** C
      PURE FUNCTION getOrientation(O) result(O2)
      ! Orders O to have the orientation x, y, z as the columns of O2
      ! @ input O: Un-ordered reference configuration
      ! @ returns O2: Ordered reference configuration such that the column of
      !               O2 [x, y, z] are a right-handed coordinate system.
      REAL(8), intent(in) :: O(3, 3)
      REAL(8)             :: O2(3, 3), xo(3), yo(3), zo(3), cross_vec(3)
      REAL(8)             :: test
      ! z is always the first column, assume that y is the second column
      zo = O(:, 1)
      yo = O(:, 2)
      xo = O(:, 3)
      ! Test to see if right-handed system
      cross_vec(1) = yo(2) * zo(3) - yo(3) * zo(2)
      cross_vec(2) = -(yo(1) * zo(3) - yo(3) * zo(1))
      cross_vec(3) = yo(1) * zo(2) - yo(2) * zo(1)
      test = DOT_PRODUCT(xo, cross_vec)
      IF (test .GT. 0) THEN
        ! Correct assumption, right hand system
        O2(:, 1) = xo
        O2(:, 2) = yo
        O2(:, 3) = zo
      ELSE
        O2(:, 1) = yo
        O2(:, 2) = xo
        O2(:, 3) = zo
      END IF
      END FUNCTION
C *************************************************************************** C
C     Compute reference configuration
C *************************************************************************** C
      ! Returns the right-hand ordered reference configuration of ref_pts.
      ! @ input ref_pts: x, y, z coordinates of each point
      ! @ returns O: Orientation of the reference configuration that is a 
      !              valid right-handed coordinate system.
      FUNCTION refConfig(ref_pts) result(O)
      !
      REAL(8), intent(in) :: ref_pts(:, :)
      REAL(8)             :: O(3, 3)
      INTEGER :: NE, NSELECT, mat_shape(2)
      PARAMETER(NE = 3, NSELECT = 3)
      INTEGER :: LDA, LDZ
      PARAMETER(LDA = NE, LDZ = NE)
      INTEGER :: INFO, LWORK, LIWORK, IL, IU, M
      DOUBLE PRECISION :: ABSTOL, VL, VU
      INTEGER :: LWMAX
      PARAMETER(LWMAX = 1000)
      INTEGER :: ISUPPZ(NE), IWORK(LWMAX)
      DOUBLE PRECISION :: AE(LDA, NE), W(NE), Z(LDZ, NSELECT), 
     1 WORK(LWMAX)
      EXTERNAL DSYEVR
      !
      mat_shape = SHAPE(ref_pts)
      AE(:, :) = ZERO
      DO i = 1, mat_shape(2)
        AE = AE 
     1  + SPREAD(ref_pts(:, i), dim=2, ncopies=3) *
     2    SPREAD(ref_pts(:, i), dim=1, ncopies=3)
      END DO
      !
      ABSTOL = -1.0
      IL = 1
      IU = NSELECT
      ! LWORK and LIWORK are set based on the query to the optimal 
      ! workspace during testing
      LWORK = 104
      LIWORK = 40
      CALL DSYEVR('Vectors', 'A', 'Upper', NE, AE, LDA, VL, VU, IL,
     1            IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     2            LIWORK, INFO) 
      ! The eigenvalues are returned in asscending order, we want the min
      ! eigenvalue to be the z-component, and the others to have a positive
      ! cross-product
      O = getOrientation(Z)
      END FUNCTION
C *************************************************************************** C
C     Calculate the warping function for all nodes
C *************************************************************************** C
      PURE FUNCTION warp_fun(C) result(psi)
      ! Returns the warping function evaluated at each node.
      ! @ input C: Location of nodes in reference config. relative to centroid.
      ! @ returns psi: Warping function each node.
      ! Input and output
      REAL(8), intent(in)   :: C(:, :)
      REAL(8), ALLOCATABLE  :: psi(:)
      ! Internal
      REAL(8)             :: s, s1, s2, TWO, FOUR, TOL
      PARAMETER              (TWO = 2.D0, FOUR = 4.D0, TOL = 1.D-8)
      INTEGER             :: i, sz(2)
      ! Allocate array
      sz = SHAPE(C)
      ALLOCATE(psi(sz(2)))
      ! Calculate the width and depth of the section
      s1 = TWO * MAXVAL(C(1, :))
      s2 = TWO * MAXVAL(C(2, :))
      DO i = 1, sz(2)
        IF (ABS(C(2, i) - s2 / TWO ) .LT. TOL) THEN
          ! Node is on the top flange
          s = C(1, i) + s1 / TWO
          psi(i) = s2 * (s1 - TWO * s) / FOUR
        ELSE IF (ABS(C(2, i) + s2 / TWO ) .LT. TOL) THEN
          ! Node is on the bottom flange
          s = C(1, i) + s1 / TWO
          psi(i) = -s2 * (s1 - TWO * s) / FOUR
        ELSE
          ! Node is on the web line
          psi(i) = 0.D0
        END IF
      END DO
      END FUNCTION
C *************************************************************************** C
C     Calculate the warping amplitude
C *************************************************************************** C
      PURE FUNCTION warp_amp(X_def, X_ref, u_c, t, R, psi) result(w)
      ! Returns the warping amplitude.
      ! @ input X_def: x of each node in the deformed configuration
      ! @ input X_ref: x of each node in the reference configuration
      ! @ input u_c: Translation of the centroid
      ! @ input t: Orientation of the normal to the cross-section
      ! @ input R: Optimal rotation matrix
      ! @ input psi: Warping function at each node
      ! @ returns w: Warping amplitude
      ! Input and output
      REAL(8), intent(in) ::  X_def(:, :), X_ref(:, :), u_c(3), t(3)
      REAL(8), intent(in) ::  R(3, 3), psi(:)
      REAL(8)             ::  w
      ! Internal
      INTEGER             :: i, non_zero_nodes, sz(2)
      ! Function start
      sz = SHAPE(X_def)
      w = 0.D0
      non_zero_nodes = 0
      DO i = 1, sz(2)
        IF (psi(i) .NE. 0.D0) THEN
          non_zero_nodes = non_zero_nodes + 1
          w = w + 1.D0 / psi(i) 
     1       *DOT_PRODUCT(t, X_def(:, i) - MATMUL(R, X_ref(:, i) + u_c))
        END IF
      END DO
      w = w / non_zero_nodes
      END FUNCTION
C *************************************************************************** C
C     Compute the linearized warping vector
C *************************************************************************** C
      PURE FUNCTION lin_warp(X_def, t0, psi, R, G) 
     1              result(w_lin)
      ! Returns the linearized warping vector.
      ! @ input x1: x at s = 0 in the deformed configuration
      ! @ input x1_ref: x at s = 0 in the reference configuration
      ! @ input u_c: Translation of the centroid
      ! @ input t0: Orientation of the normal to the cross-section in the 
      !            reference configuration
      ! @ input psi: Warping function at s = 0
      ! @ input R: Optimal rotation matrix
      ! @ input G: Instantaneous rotation matrix
      ! @ returns w_lin: Linearization of warping amplitude w.r.t X
      ! Input and output
      REAL(8), intent(in)   ::  X_def(:, :), t0(3), psi(:), R(3, 3) 
      REAL(8), intent(in)   ::  G(:, :)
      REAL(8), ALLOCATABLE  ::  w_lin(:)
      ! Internal
      INTEGER               ::  sz(2), i, j, num_nodes, num_psi_non_zero
      REAL(8)               ::  dRdX_rs(3, 3), t(3)
      REAL(8), ALLOCATABLE  ::  dRdX(:, :)
      REAL(8)               ::  ZERO, ONE, TWO
      PARAMETER                 (ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0)
      ! Allocate arrays
      ! Number of columns in G should be 3 * num_nodes
      sz = SHAPE(G)
      num_nodes = sz(2) / 3
      ALLOCATE(w_lin(sz(2)))
      ALLOCATE(dRdX(9, sz(2)))
      ! Calculate the derivative of R
      dRdX = MATMUL(-skewMat(R), G)
      ! Assemble the linearized warping vector
      t = MATMUL(R, t0)
      w_lin(:) = ZERO
      DO i = 1, sz(2)
        dRdX_rs = vec2mat_9(dRdX(:, i))
        DO j = 1, num_nodes
          IF (psi(j) .NE. ZERO) THEN
            w_lin(i) = w_lin(i) 
     1      + DOT_PRODUCT(t0, MATMUL(dRdX_rs, X_def(:, j))) / psi(j)
          END IF
        END DO
      END DO
      num_psi_non_zero = 0
      DO i = 1, num_nodes
        IF (psi(i) .NE. ZERO) THEN
          num_psi_non_zero = num_psi_non_zero + 1
          w_lin(3*i-2:3*i) = w_lin(3*i-2:3*i) + t(:) / psi(i)
        END IF
      END DO
      w_lin = w_lin / num_psi_non_zero
      END FUNCTION
C *************************************************************************** C
C     Extract rotation vector from quaternion
C *************************************************************************** C
      PURE FUNCTION extractRotation(q) result(phi)
      ! Returns the rotation vector extracted from the quaternion
      ! @ input q: Rotation quaternion
      ! @ returns phi: Euler angle representation of q
      ! From Abaqus Theory Guide 1.3.1
      ! Input and output
      REAL(8), intent(in) ::  q(4)
      REAL(8)             ::  phi(3)
      ! Internal
      REAL(8)             :: q_mag, rot
      REAL(8)             :: small
      PARAMETER              (small = 1.D-14)
      ! Function start
      q_mag = SQRT(q(2) ** 2 + q(3) ** 2 + q(4) ** 2)
      IF (q_mag .GT. small) THEN
        rot = 2.D0 * ATAN2(q_mag, q(1))
        phi = rot / q_mag * q(2:4)
      ELSE
        phi(:) = 0.D0
      END IF
      END FUNCTION
C *************************************************************************** C
      END  ! END SUBROUTINE