!234567
      SUBROUTINE get_hash(larr, darr, &
        pid, indarr, hash, hash_next)

      USE omp_lib

      IMPLICIT NONE
      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) darr(20)

      INTEGER(KIND=8) pid(larr(1))
      INTEGER(KIND=4) indarr(larr(11),2)
      INTEGER(KIND=4) hash(larr(2),larr(11)),hash_next(larr(2),larr(11))

!!!!!! Local Variables
      INTEGER(KIND=4) i, j, k, l, m
      INTEGER(KIND=4) n_thread
      INTEGER(KIND=4) nn, dn

      INTEGER(KIND=4) hash_last(larr(2),larr(11))
      INTEGER(KIND=4) ind, i0
      REAL(KIND=8) noptcl
      !REAL(KIND=8) share
      noptcl = -9223372036854775800
      nn         = larr(1)
      dn         = larr(2)
      n_thread   = larr(11)

      hash_last = -1
      !!-----
      !! MAKE HASH TABLE
      !!-----
      !$OMP PARALLEL DO default(shared) schedule(static) &
      !$OMP & private(j, ind, i0)
      DO i=1, n_thread
        DO j=indarr(i,1)+1,indarr(i,2)+1
          IF(DBLE(pid(j)) .LT. noptcl) CYCLE

          ind = MOD(ABS(pid(j)),dn) + 1
  
          IF(ind .LE. 0) ind = 1
          IF(hash(ind,i) .LT. 0) THEN
            hash(ind,i) = j
          ELSE
            i0 = hash(ind,i)
            IF(hash_last(i0,i) .LE. 0) THEN
              hash_next(i0,i) = j
              hash_last(i0,i) = j
            ELSE
              hash_next(hash_last(i0,i),i) = j
              hash_last(i0,i) = j
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE
