!234567 
      Subroutine prop_time(p_age, p_sf, p_gyr, tmp_c, tmp_s, tmp_r, tmp_g &
                , larr, darr)

      use omp_lib
      implicit none

      integer(kind=4) larr(20)
      real(kind=8) darr(20)
        
      real(kind=8) p_age(larr(1)), p_sf(larr(1)), p_gyr(larr(1))
      real(kind=8) tmp_c(10000), tmp_s(10000)
      real(kind=8) tmp_r(10000), tmp_g(10000)

      !!!!!!
      integer(kind=4) i, j, k, l
      integer(kind=4) n_part, n_thread, nn
      real(kind=8) tx(3), ty(3), aa(3), bb(3), cc(3), dum, dum0

      n_part    = larr(1)
      n_thread  = larr(2)
      nn        = 3
      call omp_set_num_threads(n_thread)

      Do j=1, 10000
        If(tmp_r(j) .gt. darr(1)) Exit
      Enddo

      tx(1) = tmp_r(j-1); tx(2) = tmp_r(j); tx(3) = tmp_r(j+1)
      ty(1) = tmp_g(j-1); ty(2) = tmp_g(j); ty(3) = tmp_g(j+1)
      call spline(tx, ty, aa, bb, cc, nn)
      dum0 = ty(2) + aa(2)*(darr(1) - tx(2)) + &
              bb(2)*(darr(1) - tx(2))**2 + &
              cc(2)*(darr(1) - tx(2))**3

      !$OMP PARALLEL DO default(shared) private(j, tx, ty, aa, bb, cc, dum) schedule(dynamic)
      Do i=1, n_part
        If(p_age(i) .gt. -1.0d7) then

          Do j=1, 10000
            If(tmp_c(j) .gt. p_age(i)) Exit
          Enddo
          tx(1) = tmp_c(j-1); tx(2) = tmp_c(j); tx(3) = tmp_c(j+1)
          ty(1) = tmp_s(j-1); ty(2) = tmp_s(j); ty(3) = tmp_s(j+1)
          call spline(tx, ty, aa, bb, cc, nn)
          p_sf(i) = ty(2) + aa(2)*(p_age(i) - tx(2)) + &
                  bb(2)*(p_age(i) - tx(2))**2 + &
                  cc(2)*(p_age(i) - tx(2))**3

          dum = 1./p_sf(i) - 1.
          Do j=1, 10000
            If(tmp_r(j) .gt. dum) Exit
          Enddo
          tx(1) = tmp_r(j-1); tx(2) = tmp_r(j); tx(3) = tmp_r(j+1)
          ty(1) = tmp_g(j-1); ty(2) = tmp_g(j); ty(3) = tmp_g(j+1)
          call spline(tx, ty, aa, bb, cc, nn)
          p_gyr(i) = ty(2) + aa(2)*(dum - tx(2)) + &
                  bb(2)*(dum - tx(2))**2 + &
                  cc(2)*(dum - tx(2))**3
          p_gyr(i) = p_gyr(i) - dum0
        Endif
      Enddo
      !$OMP END PARALLEL DO

      Return
      End subroutine

      !!!!!
      !! Spline
      !!!!!
      Subroutine spline (x, y, b, c, d, n)
      !====================================================================
      !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
      !  for cubic spline interpolation
      !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
      !  for  x(i) <= x <= x(i+1)
      !  Alex G: January 2010
      !--------------------------------------------------------------------
      !  input..
      !  x = the arrays of data abscissas (in strictly increasing order)
      !  y = the arrays of data ordinates
      !  n = size of the arrays xi() and yi() (n>=2)
      !  output..
      !  b, c, d  = arrays of spline coefficients
      !  comments ...
      !  spline.f90 program is based on fortran version of program spline.f
      !  the accompanying function fspline can be used for interpolation
      !======================================================================
      implicit none
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      integer i, j, gap
      double precision h
      
      gap = n-1
      ! check input
      if ( n < 2 ) return
      if ( n < 3 ) then
        b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
        c(1) = 0.
        d(1) = 0.
        b(2) = b(1)
        c(2) = 0.
        d(2) = 0.
        return
      end if
      !
      ! step 1: preparation
      !
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, gap
        d(i) = x(i+1) - x(i)
        b(i) = 2.0*(d(i-1) + d(i))
        c(i+1) = (y(i+1) - y(i))/d(i)
        c(i) = c(i+1) - c(i)
      end do
      !
      ! step 2: end conditions 
      !
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if(n /= 3) then
        c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)**2/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      end if
      !
      ! step 3: forward elimination 
      !
      do i = 2, n
        h = d(i-1)/b(i-1)
        b(i) = b(i) - h*d(i-1)
        c(i) = c(i) - h*c(i-1)
      end do
      !
      ! step 4: back substitution
      !
      c(n) = c(n)/b(n)
      do j = 1, gap
        i = n-j
        c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
      !
      ! step 5: compute spline coefficients
      !
      b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
      do i = 1, gap
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
        d(i) = (c(i+1) - c(i))/d(i)
        c(i) = 3.*c(i)
      end do
      c(n) = 3.0*c(n)
      d(n) = d(n-1)
      return
      end subroutine spline
