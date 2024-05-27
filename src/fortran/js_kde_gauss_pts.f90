!234567
      SUBROUTINE js_kde_gauss_pts(&
                      xx, yy, zz, xr, yr, grid, ptcl, larr)

      USE omp_lib
      IMPLICIT NONE

      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) xx(larr(1)), yy(larr(1)), zz(larr(1))
      REAL(KIND=8) xr(2), yr(2)
      REAL(KIND=8) grid(larr(3), larr(3))
      REAL(KIND=8) ptcl(larr(1))

      !LOCAL VARIABLES
      INTEGER(KIND=4) i, j, k, n_ptcl, n_thread, n_pix, ind(2)
      REAL(KIND=8) dx, dy

      n_ptcl    = larr(1)
      n_thread  = larr(2)
      n_pix     = larr(3)

      dx = (xr(2) - xr(1)) / n_pix
      dy = (yr(2) - yr(1)) / n_pix

      CALL omp_set_num_threads(n_thread)
      !$OMP PARALLEL DO default(shared) schedule(static), private(ind)
      DO i=1, n_ptcl
        ind(1) = INT((xx(i) - xr(1))/dx) + 1
        ind(2) = INT((yy(i) - yr(1))/dy) + 1

        ptcl(i) = grid(ind(1),ind(2))
      ENDDO
      !$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE
