!234567
      SUBROUTINE get_merit (larr, darr, pid0, pid, w0, w1)

      USE omp_lib

      IMPLICIT NONE
      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) darr(20)

      INTEGER(KIND=8) pid0(larr(1)), pid(larr(2))
      REAL(KIND=8) w0(larr(1)), w1(larr(2))

!!!!!! Local Variables
      INTEGER(KIND=4) i, j, k, l, m
      INTEGER(KIND=4) n_thread, n_raw, n_match

      INTEGER(KIND=8) hash(larr(1)), hash_next(larr(1))
      INTEGER(KIND=4) ind, i0
      REAL(KIND=8) noptcl
      REAL(KIND=8) share
      REAL(KIND=8) weight1, weight2
      noptcl = -9223372036854775800
      n_raw     = larr(1)
      n_match   = larr(2)
      n_thread  = larr(3)

      weight1 = 0.
      weight2 = 0.
      hash = -1
      hash_next = -1
      !!-----
      !! MAKE HASH TABLE
      !!-----
      DO i=1, n_raw
        IF(DBLE(pid0(i)).LT.noptcl) CYCLE
        ind = MOD(ABS(pid0(i)),n_raw) + 1
        IF(ind.LE.0) ind = 1
        IF(hash(ind) .LT. 0) THEN
          hash(ind) = i
        ELSE
          i0 = hash(ind)
          DO WHILE (1 .EQ. 1)
            IF(hash_next(i0) .LE. 0) THEN
              hash_next(i0) = i
              EXIT
            ELSE
              i0 = hash_next(i0)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      !!-----
      !! MATCHING
      !!-----

      CALL OMP_SET_NUM_THREADS(n_thread)
      share = 0.
      !$OMP PARALLEL DO default(shared) schedule(static) &
      !$OMP & private(ind, i0) reduction(+:share, weight1, weight2)
      DO i=1, n_match
        IF(DBLE(pid(i)) .LE. noptcl) CYCLE
        ind = MOD(ABS(pid(i)), n_raw) + 1
        IF(ind.LE.0) ind = 1
        i0 = hash(ind)
        IF(i0 .LE. 0) CYCLE
        DO WHILE(1 .EQ. 1)
          IF(pid(i) .EQ. pid0(i0)) THEN
            share = share + 1
            weight1 = weight1 + w0(i0)
            weight2 = weight2 + w1(i)
            EXIT
          ELSE
            IF(hash_next(i0) .GT. 0) THEN
              i0 = hash_next(i0)
            ELSE
              EXIT
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      IF( larr(4) .EQ. 1 ) darr(1) = share * share / DBLE(n_raw) / DBLE(n_match) 
      IF( larr(4) .EQ. 2 ) darr(1) = share / DBLE(n_raw)
      IF( larr(4) .EQ. 3 ) darr(1) = share * share / DBLE(n_raw) / DBLE(n_match) * weight1 * weight2
      END SUBROUTINE
