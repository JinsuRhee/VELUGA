!234567
      SUBROUTINE js_kde_gauss_binning(&
                      xx, yy, zz, xr, yr, grid, larr)

      USE omp_lib
      IMPLICIT NONE

      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) xx(larr(1)), yy(larr(1)), zz(larr(1))
      REAL(KIND=8) xr(2), yr(2)
      REAL(KIND=8) grid(larr(3)*2, larr(3)*2)

      !LOCAL VARIABLES
      INTEGER(KIND=4) i, j, k, n_ptcl, n_thread, n_pix, n_pix2
      INTEGER(KIND=4) bintype, histtype
      INTEGER(KIND=4) ind(larr(1),2)
      REAL(KIND=8) xr2(2), yr2(2), dx, dy, indval(larr(1),2)
      REAL(KIND=8) bindata(larr(1),4)

      n_ptcl    = larr(1)
      n_thread  = larr(2)
      n_pix     = larr(3)
      bintype   = larr(4)
      histtype  = larr(5)

      n_pix2    = n_pix*2

      xr2(1) = xr(1) - (xr(2) - xr(1))*0.5
      xr2(2) = xr(2) + (xr(2) - xr(1))*0.5
      yr2(1) = yr(1) - (yr(2) - yr(1))*0.5
      yr2(2) = yr(2) + (yr(2) - yr(1))*0.5

      dx = (xr(2) - xr(1)) / n_pix
      dy = (yr(2) - yr(1)) / n_pix

      CALL omp_set_num_threads(n_thread)
      !$OMP PARALLEL DO default(shared) schedule(static,100)
      DO i=1, n_ptcl
        ind(i,1) = INT((xx(i) - xr2(1))/dx) + 1
        ind(i,2) = INT((yy(i) - yr2(1))/dy) + 1

        indval(i,1) = xr2(1) + dx*(ind(i,1) - 0.5)
        indval(i,2) = yr2(1) + dy*(ind(i,2) - 0.5)

        IF(indval(i,1) .GT. xx(i)) THEN
          ind(i,1) = ind(i,1) - 1
          indval(i,1) = indval(i,1) - dx
        ENDIF

        IF(indval(i,2) .GT. yy(i)) THEN
          ind(i,2) = ind(i,2) - 1
          indval(i,2) = indval(i,2) - dy
        ENDIF

        IF(bintype .EQ. 1) THEN         ! NGP
                bindata(i,1) = zz(i)
                bindata(i,2) = 0.
                bindata(i,3) = 0.
                bindata(i,4) = 0.
        ELSE IF(bintype .EQ. 2) THEN    ! CIC
                bindata(i,1) = zz(i)*1./4 * &
                        (dx - (xx(i) - indval(i,1))) / dx * &
                        (dy -(yy(i) - indval(i,2)))/dy
                bindata(i,2) = zz(i)*1./4 * &
                        (xx(i) - indval(i,1)) / dx * &
                        (dy - (yy(i) - indval(i,2))) / dy
                bindata(i,3) = zz(i)*1./4 * &
                        (dx - (xx(i) - indval(i,1))) / dx * &
                        (yy(i) - indval(i,2)) / dy
                bindata(i,4) = zz(i)*1./4 * (&
                        xx(i) - indval(i,1)) / dx * &
                        (yy(i) - indval(i,2)) / dy
        ENDIF
      ENDDO
      !$OMP END PARALLEL DO

      DO i=1, n_ptcl
        j = ind(i,1)
        k = ind(i,2)

        IF(j.GT.n_pix*2) j = n_pix*2
        IF(k.GT.n_pix*2) k = n_pix*2

        IF(histtype .EQ. 1) THEN                ! Normal HISTO
                grid(j  ,k  ) = grid(j  ,k  ) + bindata(i,1)
                !grid(j+1,k  ) = grid(j+1,k  ) + bindata(i,2)
                !grid(j  ,k+1) = grid(j  ,k+1) + bindata(i,3)
                !grid(j+1,k+1) = grid(j+1,k+1) + bindata(i,4)
        ELSE IF(histtype .EQ. 2)THEN            ! Maximum
                grid(j  ,k  ) = MAX(grid(j  ,k  ), bindata(i,1))
                grid(j+1,k  ) = MAX(grid(j+1,k  ), bindata(i,2))
                grid(j  ,k+1) = MAX(grid(j  ,k+1), bindata(i,3))
                grid(j+1,k+1) = MAX(grid(j+1,k+1), bindata(i,4))
        END IF
      ENDDO

      RETURN
      END SUBROUTINE
