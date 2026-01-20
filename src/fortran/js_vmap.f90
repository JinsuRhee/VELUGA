!234567
      SUBROUTINE js_vmap(larr, darr, ix, iy, vx, vy, mass, vmap)

      USE omp_lib
      IMPLICIT NONE
      REAL(KIND=8) darr(20)
      INTEGER(KIND=4) larr(20)

      INTEGER(KIND=4) ix(larr(1)), iy(larr(1))
      REAL(KIND=8) vx(larr(2)), vy(larr(2)), mass(larr(2))
      REAL(KIND=8) vmap(larr(3),larr(3),3)


      !!!!!
      !! LOCAL
      !!!!!
      INTEGER(KIND=4) i, j, k
      INTEGER(KIND=4) n_thread, n_i, n_v, n_pix

      n_i = larr(1)
      n_v = larr(2)
      n_pix = larr(3)
      n_thread = larr(4)

      CALL OMP_SET_NUM_THREADS(n_thread)

      ix = ix + 1
      iy = iy + 1

      !$OMP PARALLEL DO default(shared) &
      !$OMP & schedule(static) &
      !$OMP & reduction(+:vmap)
      DO i=1, n_v
        IF(ix(i) .LT. 1 .OR. ix(i) .GT. n_pix) CYCLE
        IF(iy(i) .LT. 1 .OR. iy(i) .GT. n_pix) CYCLE

        vmap(ix(i),iy(i),1) = vmap(ix(i),iy(i),1) + vx(i) * mass(i)
        vmap(ix(i),iy(i),2) = vmap(ix(i),iy(i),2) + vy(i) * mass(i)
        vmap(ix(i),iy(i),3) = vmap(ix(i),iy(i),3) + mass(i)
      ENDDO
      !$OMP END PARALLEL DO

      END SUBROUTINE

