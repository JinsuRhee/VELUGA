!234567
      SUBROUTINE get_merit2(larr, darr, &
        pid_g, gid_g, pid_s, gid_s, &
        hash, hash_next, &
        npart_g, npart_s, &
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
      !REAL(KIND=8) merit(larr(3), larr(4))

      INTEGER(KIND=4) m_id(larr(7))
      REAL(KIND=8) m_merit(larr(7))

!!!!!! Local Variables
      INTEGER(KIND=4) i, j, k, l, m
      INTEGER(KIND=4) n_thread, n_pg, n_ps
      INTEGER(KIND=4) n_g, n_s, n_tree
      INTEGER(KIND=4) n_dn, n_hashth
      INTEGER(KIND=4) idtoind(larr(3)), offomp

      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: merit
      INTEGER(KIND=8) tbl_g(larr(3)), tbl_s(larr(4))
      INTEGER(KIND=4) nn_g, nn_s

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
      offomp     = larr(12)

      !!-----
      !! MAKE INDEX TABLE
      !!-----
      tbl_g = 0
      tbl_s = 0
      nn_g = 1
      nn_s = 1
!PRINT *, 'AA'
      DO i=1, larr(3)
        IF(npart_g(i) .GT. 0) THEN
                tbl_g(i) = nn_g
                nn_g = nn_g + 1
        ENDIF
      ENDDO

      DO i=1, larr(4)
        IF(npart_s(i) .GT. 0) THEN
                tbl_s(i) = nn_s
                nn_s = nn_s + 1
        ENDIF
      ENDDO
      nn_g = nn_g - 1
      nn_s = nn_s - 1

      ALLOCATE(merit(1:nn_g, 1:nn_s))
!PRINT *, 'BB' 
      !!-----
      !! MATCHING
      !!-----
      IF(offomp .LT. 0) THEN
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
                !merit(gid_g(i)+1, gid_s(i0)+1) = merit(gid_g(i)+1,gid_s(i0)+1) + 1
                merit(tbl_g(gid_g(i)+1), tbl_s(gid_s(i0)+1)) = merit(tbl_g(gid_g(i)+1), tbl_s(gid_s(i0)+1)) + 1
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
      ELSE
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
                !merit(gid_g(i)+1, gid_s(i0)+1) = merit(gid_g(i)+1,gid_s(i0)+1) + 1
                !PRINT *, i, j, gid_g(i), gid_s(i0)
                !PRINT *,        tbl_g(gid_g(i))
                !PRINT *,        tbl_s(gid_s(i0))
                merit(tbl_g(gid_g(i)+1), tbl_s(gid_s(i0)+1)) = merit(tbl_g(gid_g(i)+1), tbl_s(gid_s(i0)+1)) + 1
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
      ENDIF
      merit     = merit * merit
!PRINT *, 'CC'
      !$OMP PARALLEL DO default(shared) schedule(static)
      DO i=1, n_g
        IF(npart_g(i) .GT. 0) merit(tbl_g(i),:) = merit(tbl_g(i),:) / npart_g(i)
      ENDDO
      !$OMP END PARALLEL DO
!PRINT *, 'DD'
      !$OMP PARALLEL DO default(shared) schedule(static)
      DO i=1, n_s
        IF(npart_s(i) .GT. 0) merit(:,tbl_s(i)) = merit(:,tbl_s(i)) / npart_s(i)
      ENDDO
      !$OMP END PARALLEL DO
      
      ind = 0
      DO i=1, n_g
        IF(npart_g(i) .GT. 0) THEN
          ind = ind + 1
          idtoind(i) = ind
        ENDIF
      ENDDO
!PRINT *, 'EE'
      !$OMP PARALLEL DO default(shared) schedule(static) &
      !$OMP & private(dum, dumid, j)
      DO i=1, n_g
        IF(npart_g(i) .LE. 0) CYCLE

        dum = -1.0
        dumid = -1

        DO j=1, n_s
          IF(npart_s(j) .LE. 0) CYCLE
          IF(merit(tbl_g(i),tbl_s(j)) .GT. dum .AND. merit(tbl_g(i),tbl_s(j)) .GT. 0) THEN
            dum = merit(tbl_g(i),tbl_s(j))
            dumid = j-1
          ENDIF
        ENDDO

        IF(dumid .GT. 0) THEN
          m_id(idtoind(i)) = dumid
          m_merit(idtoind(i)) = dum
        ENDIF

      ENDDO
      !$OMP END PARALLEL DO
!PRINT *, 'FF'
      DEALLOCATE(merit)
      END SUBROUTINE
