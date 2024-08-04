!234567
      SUBROUTINE js_indmatch(larr, darr, x, y, x_m, y_m)

      USE omp_lib
      IMPLICIT NONE
      REAL(KIND=8) darr(20)
      INTEGER(KIND=4) larr(20)

      INTEGER(KIND=8) x(larr(2)), y(larr(3))
      INTEGER(KIND=8) x_m(larr(2)), y_m(larr(3))

      !!-----
      !! LOCAL
      !!-----
      INTEGER(KIND=4) i, j, k
      INTEGER(KIND=4) num_thread, nx, ny, stype

      num_thread = larr(1)
      nx = larr(2)
      ny = larr(3)
      stype = larr(11)
      CALL omp_set_num_threads(num_thread)
      IF(stype .EQ. 1) CALL js_indmatch_hash(x, y, x_m, y_m, nx, ny)

      CONTAINS
!!-----
!! FIND INDEX
!!-----
      SUBROUTINE js_indmatch_hash(x, y, x_m, y_m, nx, ny)
      INTEGER(KIND=4) nx, ny, i, j, k
      INTEGER(KIND=4) h_tbl(ny), h_tind(ny)
      INTEGER(KIND=4) h_ind(ny), h_ind2(ny), h_ind3(ny)
      INTEGER(KIND=8) x(nx), y(ny), x_m(nx), y_m(ny)

      !! MAKE HASH TABLE WITH Y
      h_ind = 0
      !! k
      !! reduced sum
      DO i=1, ny
        k = MOD(y(i), ny) + 1
        h_ind(k) = h_ind(k) + 1
      ENDDO


      h_ind2(1) = 1
      DO i=2, ny
        h_ind2(i) = h_ind(i-1) + h_ind2(i-1)
      ENDDO

      h_ind3 = h_ind2
      DO i=1, ny
        k = MOD(y(i),ny) + 1
        h_tbl(h_ind3(k)) = y(i)
        h_tind(h_ind3(k)) = i
        h_ind3(k) = h_ind3(k) + 1
      ENDDO

      !! SEARCH BY HASH
      DO i=1, nx
        k = MOD(x(i), ny) + 1

        IF(h_ind(k) .EQ. 0) CYCLE

        DO j=h_ind2(k), h_ind2(k) + h_ind(k)-1
          IF(h_tbl(j) .EQ. x(i)) THEN
            x_m(i) = h_tind(j) - 1
            y_m(h_tind(j)) = i - 1
            EXIT
          ENDIF
        ENDDO

      ENDDO
      END SUBROUTINE js_indmatch_hash

      END SUBROUTINE js_indmatch


