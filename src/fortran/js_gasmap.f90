!234567
      SUBROUTINE js_gasmap(larr, darr, xx, yy, zz, bw, xr, yr, map, amrtype)

      USE omp_lib
      IMPLICIT NONE
      REAL(KIND=8) darr(20)
      INTEGER(KIND=4) larr(20)

      REAL(KIND=8) xx(larr(1)), yy(larr(1)), zz(larr(1),1+larr(2)+1)
      REAL(KIND=8) bw(2), xr(2), yr(2)
      REAL(KIND=8) map(larr(3), larr(3), 2)
      INTEGER(KIND=4) amrtype(larr(2))

      !!!!!
      !! LOCAL
      !!!!!
      INTEGER(KIND=4) i, j, k, l
      INTEGER(KIND=4) n_thread, n_cell, n_pix, n_amr
      REAL(KIND=8) dx, dy, test, geometry
      INTEGER(KIND=4) nrx, nry, nx0, nx1, ny0, ny1, nx, ny
      !INTEGER(KIND=4) amrtype

      n_cell = larr(1)
      n_amr = larr(2)
      n_pix  = larr(3)
      n_thread= larr(4)
      !amrtype=larr(11)

      dx = (xr(2) - xr(1))/n_pix
      dy = (yr(2) - yr(1))/n_pix
      nrx = INT(bw(1)/dx/2.) + 1
      nry = INT(bw(2)/dy/2.) + 1

      CALL OMP_SET_NUM_THREADS(n_thread)

      !IF(amrtype .NE. 3) THEN
          !!$OMP PARALLEL DO &
          !!$OMP & shared(xx, yy, xr, yr, dx, dy, nrx, nry, bw) &
          !!$OMP & shared(amrtype, zz, n_pix) &
          !!$OMP & private(nx, ny, nx0, nx1, ny0, ny1, j, k, geometry) &
          !!$OMP & reduction(+:map)
          DO i=1, n_cell
            IF(xx(i) .LT. xr(1)) CYCLE
            IF(xx(i) .GT. xr(2)) CYCLE
            IF(yy(i) .LT. yr(1)) CYCLE
            IF(yy(i) .GT. yr(2)) CYCLE
    
            nx = INT((xx(i) - xr(1))/dx) + 1
            ny = INT((yy(i) - yr(1))/dy) + 1
    
            IF(nrx .GE. 1) THEN
                    nx0 = MAX(nx-nrx,1)
                    nx1 = MIN(nx+nrx,n_pix)
            ELSE
                    nx0 = nx
                    nx1 = nx
            ENDIF
    
            IF(nry .GE. 1) THEN
                    ny0 = MAX(ny-nry,1)
                    ny1 = MIN(ny+nry,n_pix)
            ELSE
                    ny0 = ny
                    ny1 = ny
            ENDIF
    
            DO j=nx0, nx1
            DO k=ny0, ny1
              !compute partial volume
              !geometry*bw(1) = partial volume
              CALL grid_geometry(j, k, dx, dy, xr(1), yr(1), xx(i), yy(i), bw(1), bw(2), geometry)
    

              !Store mass map
              map(j,k,1)    = map(j,k,1) + zz(i,1) * geometry * bw(1)
    
              DO l=1, n_amr
                IF(amrtype(l) .EQ. 1) THEN !! Mass-weighted
                        map(j,k,l+1) = map(j,k,l+1) + zz(i,l+1) * zz(i,1) * geometry * bw(1)
                ELSE IF (amrtype(l) .EQ. 2) THEN !! Volume-weighted
                        map(j,k,l+1) = map(j,k,l+1) + zz(i,l+1) * geometry * bw(1)
                ELSE IF (amrtype(l) .EQ. 3) THEN !! MAX
                        !map(j,k,1) = MAX(map(j,k,1), geometry * zz(i,1))
                        map(j,k,l+1) = MAX(map(j,k,l+1), zz(i,l+1))
                ELSE IF (amrtype(l) .EQ. 4) THEN !! Column Density
                        map(j,k,l+1)    = map(j,k,l+1) + zz(i,1)*geometry*bw(1) * darr(1)
                        !map(j,k,1) = map(j,k,1) + zz(i,2)*geometry*(xr(2)-xr(1))
                        !map(j,k,1) = map(j,k,1) + zz(i,2) * geometry * bw(1) / (dx*dy)
                        !map(j,k,1) = map(j,k,1) + zz(i,2) * geometry * bw(1) * zz(i,2)
                        !map(j,k,1) = map(j,k,1) + bw(1)!zz(i,1) * geometry *bw(1)**3
                ELSE IF (amrtype(l) .EQ. 5) THEN !! HISTOGRAM
                        map(j,k,l+1) = map(j,k,l+1) + zz(i,l+1)*(geometry/(bw(1)*bw(2)))
                ENDIF
              ENDDO
    
            ENDDO
            ENDDO
            !  IF(amrtype .EQ. 1) THEN !! TEMPERATURE
            !    IF(zz(i,1) .GT. map(j,k)) map(j,k) = zz(i,1)
            !  ENDIF
            !  IF(amrtype .EQ. 2) THEN !! DENSITY
            !        map(j,k) = map(j,k) + zz(i,2)*bw(1)
            !  ENDIF
            !  IF(amrtype .EQ. 3) THEN !! DENSITY WEIGHTED TEMPERATURE
            !    IF(zz(i,1)*zz(i,2)*bw(1) .GT. map(j,k)) map(j,k) = zz(i,1)*zz(i,2)*bw(1)
            !        !map(j,k) = map(j,k) + zz(i,1) * zz(i,2)*bw(1)
            !  ENDIF
            !ENDDO
            !ENDDO
          ENDDO
          !!$OMP END PARALLEL DO
      !ELSE
          !!$OMP PARALLEL DO &
          !!$OMP & shared(xx, yy, xr, yr, dx, dy, nrx, nry, bw) &
          !!$OMP & shared(amrtype, zz) &
          !!$OMP & private(nx, ny, nx0, nx1, ny0, ny1, j, k, geometry) &
          !!$OMP & reduction(max:map)
          
          !!$OMP END PARALLEL DO
      !ENDIF

      END SUBROUTINE

      SUBROUTINE grid_geometry(ix, iy, dx, dy, xr0, yr0, x, y, bx, by, geometry)

      IMPLICIT NONE
      INTEGER(KIND=4) ix, iy, i, j, k
      REAL(KIND=8) dx, dy, x, y, bx, by, geometry
      REAL(KIND=8) x0, y0, x1, x2, y1, y2, delX, delY, xr0, yr0


      geometry = 0

      x0 = (ix-0.5)*dx + xr0
      y0 = (iy-0.5)*dy + yr0

      x1 = x - bx/2.
      x2 = x + bx/2.
      y1 = y - by/2.
      y2 = y + by/2.

      IF(x0 + dx/2. .LT. x1) THEN
        RETURN
      ELSE IF(x0 - dx/2. .GT. x2) THEN
        RETURN
      ELSE IF(y0 + dy/2. .LT. y1) THEN
        RETURN
      ELSE IF(y0 - dy/2. .GT. y2) THEN
        RETURN
      ELSE IF(ABS(x0-x1) .LT. dx/2.) THEN
        delX = x0 + dx/2. - x1
        IF(ABS(y0-y1) .LT. dy/2.) THEN
          delY = y0 + dy/2. - y1
        ELSE IF(ABS(y0-y2) .LT. dy/2.) THEN
          delY = y2 - (y0 - dy/2.)
        ELSE
          delY = dy
        ENDIF
        geometry = delX * delY
        RETURN
      ELSE IF(ABS(x0-x2) .LT. dx/2.) THEN
        delX = x2 - (x0 - dx/2.)
        IF(ABS(y0-y1) .LT. dy/2.) THEN
          delY = y0 + dy/2. - y1
        ELSE IF(ABS(y0-y2) .LT. dy/2.) THEN
          delY = y2 - (y0 - dy/2.)
        ELSE
          delY = dy
        ENDIF
        geometry = delX * delY
        RETURN
      ELSE IF(ABS(y0-y1) .LT. dy/2.) THEN
        delY = y0 + dy/2. - y1
        delX = dx
        geometry = delX * delY
        RETURN
      ELSE IF(ABS(y0-y2) .LT. dy/2.) THEN
        delY = y2 - (y0 - dy/2.)
        delX = dx
        geometry = delX * delY
        RETURN
      ELSE
        geometry = dx * dY
        RETURN
      ENDIF
      END SUBROUTINE

