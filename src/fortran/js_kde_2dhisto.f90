!234567
      Subroutine js_kde_2dhisto(xx, yy, zz, grid, ptcl, xr, yr, larr)

      Use omp_lib
      Implicit none
      INTEGER(KIND=4) larr(20)

      REAL(KIND=8) xx(larr(1)), yy(larr(1)), zz(larr(1))
      REAL(KIND=8) xr(2), yr(2)
      REAL(KIND=4) grid(larr(3), larr(3)), ptcl(larr(1))

      !!!!!
      !!!!!

      INTEGER(KIND=4) i, j, k, num_thread, n_part, ind(larr(1),2), n_pix
      REAL(KIND=8) dx, dy

      CALL omp_set_num_threads(num_thread)

      !!!!!
      !! Settings
      !!!!!

      n_part = larr(1)
      num_thread = larr(2)
      n_pix = larr(3)

      dx = (xr(2) - xr(1)) / n_pix
      dy = (yr(2) - yr(1)) / n_pix

      !$OMP PARALLEL DO default(shared) schedule(static)
      DO i=1, n_part
        ind(i,1) = INT((xx(i) - xr(1))/dx) + 1
        ind(i,2) = INT((yy(i) - yr(1))/dy) + 1
      ENDDO
      !$OMP END PARALLEL DO

      DO i=1, n_part
        grid(ind(i,1),ind(i,2)) = grid(ind(i,1),ind(i,2)) + zz(i)
      ENDDO

      !$OMP PARALLEL DO default(shared) schedule(static)
      DO i=1, n_part
        ptcl(i) = grid(ind(i,1),ind(i,2))
      ENDDO
      !$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE
