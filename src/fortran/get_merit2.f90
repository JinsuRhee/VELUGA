!234567
      SUBROUTINE get_merit2(larr, darr, &
        pid_g, gid_g, pid_s, gid_s, &
        hash, hash_next, &
        npart_g, npart_s, merit, &
        m_id, m_merit)

      USE omp_lib

      IMPLICIT NONE
      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) darr(20)

      INTEGER(KIND=8) pid_g(larr(1)), pid_s(larr(2))
      INTEGER(KIND=4) gid_g(larr(1)), gid_s(larr(2))
      INTEGER(KIND=4) npart_g(larr(3)), npart_s(larr(4))
      INTEGER(KIND=4) hash(larr(5),larr(6))
      INTEGER(KIND=4) hash_next(larr(5),larr(6))
      REAL(KIND=8) merit(larr(3), larr(4))

      INTEGER(KIND=4) m_id(larr(7))
      REAL(KIND=8) m_merit(larr(7))

!!!!!! Local Variables
      INTEGER(KIND=4) i, j, k, l, m
      INTEGER(KIND=4) n_thread, n_pg, n_ps
      INTEGER(KIND=4) n_g, n_s, n_tree
      INTEGER(KIND=4) n_dn, n_hashth
      INTEGER(KIND=4) idtoind(larr(3))

      INTEGER(KIND=4) ind, i0, check, dumid
      REAL(KIND=8) noptcl, dum
      !REAL(KIND=8) share
      noptcl = -9223372036854775800
      n_pg       = larr(1)
      n_ps       = larr(2)
      n_g        = larr(3)
      n_s        = larr(4)
      n_dn       = larr(5)
      n_hashth   = larr(6)
      n_tree     = larr(7)
      n_thread   = larr(11)
      !!-----
      !! MATCHING
      !!-----

      CALL OMP_SET_NUM_THREADS(n_thread)
      !$OMP PARALLEL DO default(shared) schedule(static) &
      !$OMP & private(ind, i0, j, check) reduction(+:merit)
      DO i=1, n_pg
        IF(DBLE(pid_g(i)) .LE. noptcl) CYCLE
        ind = MOD(ABS(pid_g(i)), n_dn) + 1
        IF(ind.LE.0) ind = 1

        !!
        check = -1
        DO j=1, n_hashth
          i0 = hash(ind,j)
          IF(i0 .LE. 0) CYCLE
          DO WHILE(1 .EQ. 1)
            IF(pid_g(i) .EQ. pid_s(i0)) THEN
              merit(gid_g(i), gid_s(i0)) = merit(gid_g(i), gid_s(i0)) + 1
              darr(1) = darr(1) + 1
              check = 1
              EXIT
            ELSE
              IF(hash_next(i0,j) .GT. 0) THEN
                i0 = hash_next(i0,j)
              ELSE
                EXIT
              ENDIF
            ENDIF
          ENDDO

          IF(check .GT. 0) EXIT
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      merit     = merit * merit

      !$OMP PARALLEL DO default(shared) schedule(static)
      DO i=1, n_g
        IF(npart_g(i) .GT. 0) merit(i,:) = merit(i,:) / npart_g(i)
      ENDDO
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO default(shared) schedule(static)
      DO i=1, n_s
        IF(npart_s(i) .GT. 0) merit(:,i) = merit(:,i) / npart_s(i)
      ENDDO
      !$OMP END PARALLEL DO

      ind = 0
      DO i=1, n_g
        IF(npart_g(i) .GT. 0) THEN
          ind = ind + 1
          idtoind(i) = ind
        ENDIF
      ENDDO

      !$OMP PARALLEL DO default(shared) schedule(static) &
      !$OMP & private(dum, dumid, j)
      DO i=1, n_g
        IF(npart_g(i) .LE. 0) CYCLE

        dum = -1.0
        dumid = -1

        DO j=1, n_s
          IF(npart_s(j) .LE. 0) CYCLE
          IF(merit(i,j) .GT. dum .AND. merit(i,j) .GT. 0) THEN
            dum = merit(i,j)
            dumid = j
          ENDIF
        ENDDO

        IF(dumid .GT. 0) THEN
          m_id(idtoind(i)) = dumid
          m_merit(idtoind(i)) = dum
        ENDIF
      ENDDO
      !$OMP END PARALLEL DO


      END SUBROUTINE
