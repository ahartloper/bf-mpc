      subroutine jacobi_eigenvalue(n,a,it_max,v,d,it_num,rot_num)
!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
!
!    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
!
!    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
!
      implicit none

      integer ( kind = 4 ) n

      real ( kind = 8 ) a(n,n)
      real ( kind = 8 ) bw(n)
      real ( kind = 8 ) c
      real ( kind = 8 ) d(n)
      real ( kind = 8 ) g
      real ( kind = 8 ) gapq
      real ( kind = 8 ) h
      integer ( kind = 4 ) i
      integer ( kind = 4 ) it_max
      integer ( kind = 4 ) it_num
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      integer ( kind = 4 ) l
      integer ( kind = 4 ) m
      integer ( kind = 4 ) p
      integer ( kind = 4 ) q
      integer ( kind = 4 ) rot_num
      real ( kind = 8 ) s
      real ( kind = 8 ) t
      real ( kind = 8 ) tau
      real ( kind = 8 ) term
      real ( kind = 8 ) termp
      real ( kind = 8 ) termq
      real ( kind = 8 ) theta
      real ( kind = 8 ) thresh
      real ( kind = 8 ) v(n,n)
      real ( kind = 8 ) w(n)
      real ( kind = 8 ) zw(n)

      do j = 1, n
        do i = 1, n
          v(i,j) = 0.0D+00
        end do
        v(j,j) = 1.0D+00
      end do

      do i = 1, n
        d(i) = a(i,i)
      end do

      bw(1:n) = d(1:n)
      zw(1:n) = 0.0D+00
      it_num = 0
      rot_num = 0

      do while ( it_num < it_max )

        it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
      thresh = 0.0D+00
      do j = 1, n
        do i = 1, j - 1
          thresh = thresh + a(i,j) ** 2
        end do
      end do

      thresh = sqrt ( thresh ) / real ( 4 * n, kind = 8 )

      if ( thresh == 0.0D+00 ) then
        exit 
      end if

      do p = 1, n
        do q = p + 1, n

          gapq = 10.0D+00 * abs ( a(p,q) )
          termp = gapq + abs ( d(p) )
          termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. 
     1        termp == abs ( d(p) ) .and. 
     2        termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) 
     1       + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then 
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h                  
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
              do j = 1, n
                g = v(j,p)
                h = v(j,q)
                v(j,p) = g - s * ( h + g * tau )
                v(j,q) = h + s * ( g - h * tau )
              end do

              rot_num = rot_num + 1

            end if

          end do
        end do

        bw(1:n) = bw(1:n) + zw(1:n)
        d(1:n) = bw(1:n)
        zw(1:n) = 0.0D+00

      end do
!
!  Restore upper triangle of input matrix.
!
      do j = 1, n
        do i = 1, j - 1
          a(i,j) = a(j,i)
        end do
      end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
        do k = 1, n - 1

          m = k

          do l = k + 1, n
            if ( d(l) < d(m) ) then
              m = l
            end if
          end do

          if ( m /= k ) then

            t    = d(m)
            d(m) = d(k)
            d(k) = t

            w(1:n)   = v(1:n,m)
            v(1:n,m) = v(1:n,k)
            v(1:n,k) = w(1:n)

          end if

        end do

        return
      end


      subroutine r8mat_diag_get_vector ( n, a, v )
!*****************************************************************************80
!
!! R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) V(N), the diagonal entries
!    of the matrix.
!
        implicit none

        integer ( kind = 4 ) n

        real ( kind = 8 ) a(n,n)
        integer ( kind = 4 ) i
        real ( kind = 8 ) v(n)

        do i = 1, n
          v(i) = a(i,i)
        end do

        return
      end


      subroutine r8mat_identity ( n, a )
!*****************************************************************************80
!
!! R8MAT_IDENTITY stores the identity matrix in an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N identity matrix.
!
      implicit none

      integer ( kind = 4 ) n

      real ( kind = 8 ) a(n,n)
      integer ( kind = 4 ) i

      a(1:n,1:n) = 0.0D+00

      do i = 1, n
        a(i,i) = 1.0D+00
      end do

      return
      end


      subroutine r8mat_is_eigen_right ( n, k, a, x, lambda, 
     1   error_frobenius )
!*****************************************************************************80
!
!! R8MAT_IS_EIGEN_RIGHT determines the error in a (right) eigensystem.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    This routine computes the Frobenius norm of
!
!      A * X - X * LAMBDA
!
!    where
!
!      A is an N by N matrix,
!      X is an N by K matrix (each of K columns is an eigenvector)
!      LAMBDA is a K by K diagonal matrix of eigenvalues.
!
!    This routine assumes that A, X and LAMBDA are all real!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) K, the number of eigenvectors.
!    K is usually 1 or N.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N,K), the K eigenvectors.
!
!    Input, real ( kind = 8 ) LAMBDA(K), the K eigenvalues.
!
!    Output, real ( kind = 8 ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
!    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
!
      implicit none

      integer ( kind = 4 ) k
      integer ( kind = 4 ) n

      real ( kind = 8 ) a(n,n)
      real ( kind = 8 ) c(n,k)
      real ( kind = 8 ) error_frobenius
      integer ( kind = 4 ) j
      real ( kind = 8 ) lambda(k)
      real ( kind = 8 ) r8mat_norm_fro
      real ( kind = 8 ) x(n,k)

      c(1:n,1:k) = matmul ( a(1:n,1:n), x(1:n,1:k) )

      do j = 1, k
        c(1:n,j) = c(1:n,j) - lambda(j) * x(1:n,j)
      end do

      error_frobenius = r8mat_norm_fro ( n, k, c )

      return
      end


      function r8mat_norm_fro ( m, n, a )
!*****************************************************************************80
!
!! R8MAT_NORM_FRO returns the Frobenius norm of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = 8 ) values.
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) A(I,J)^2 )
!
!    The matrix Frobenius-norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A*x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_FRO, the Frobenius norm of A.
!
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      real ( kind = 8 ) a(m,n)
      real ( kind = 8 ) r8mat_norm_fro

      r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

      return
      end
