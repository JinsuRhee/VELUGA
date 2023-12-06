!234567
      SUBROUTINE jsrd_part(larr, darr, fname, part_ind, &
        xp, vp, mp, ap, zp, fam, tag, domain, idvar, domlist)
        !dblvar, lonvar, idvar, domlist)

      USE omp_lib
      IMPLICIT NONE
      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) darr(20)

      CHARACTER(larr(4)) fname
      INTEGER(KIND=4) part_ind(larr(1))
      REAL(KIND=8) xp(larr(5),3), vp(larr(5),3)
      REAL(KIND=8) mp(larr(5)), ap(larr(5)), zp(larr(5))
      INTEGER(KIND=4) fam(larr(5)), tag(larr(5)), domain(larr(5))
      !REAL(KIND=8) dblvar(larr(5),9)
      !INTEGER(KIND=4) lonvar(larr(5),3)
      INTEGER(KIND=8) idvar(larr(5))
      INTEGER(KIND=4) domlist(larr(1))

!!!!! LOCAL VARIABLES
      INTEGER(KIND=4) i, j, k, i2
      INTEGER(KIND=4) n_thread, cpu0, cpu1 , uout, npart, ncpu

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tmp_dbl
      INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE :: tmp_byte
      INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: tmp_int
      INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: tmp_l64
      CHARACTER*(100) domnum, fdum
      ncpu = larr(1)
      n_thread = larr(3)

      IF (n_thread .GT. ncpu) n_thread=ncpu
      CALL OMP_SET_NUM_THREADS(n_thread)

      !domnum, fdum, uout, npart, tmp_dbl, tmp_l64, tmp_int, tmp_byte, j, k
      !$OMP PARALLEL DO default(shared) &
      !$OMP & schedule(static,10) private(domnum, fdum, uout, npart) &
      !$OMP & private(tmp_dbl, tmp_l64, tmp_int, tmp_byte, j, k, i)
      DO i2=1, ncpu
        i = domlist(i2)
        WRITE(domnum, '(I5.5)') i
        fdum = TRIM(fname)//TRIM(domnum)

        uout = OMP_GET_THREAD_NUM()
        uout = uout + 10

        OPEN(UNIT=uout, FILE=fdum, FORM='unformatted', STATUS='old')
        READ(uout); READ(uout); READ(uout) npart; READ(uout);
        READ(uout); READ(uout); READ(uout); READ(uout);

        ALLOCATE(tmp_dbl(1:npart))
        ALLOCATE(tmp_l64(1:npart))
        ALLOCATE(tmp_int(1:npart))
        ALLOCATE(tmp_byte(1:npart))

        !x, y, z
        DO j=1, 3
          READ(uout) tmp_dbl
          DO k=1, npart
            !dblvar(part_ind(i2)+k,j) = tmp_dbl(k)
            xp(part_ind(i2)+k,j) = tmp_dbl(k)
          ENDDO
        ENDDO

        !vx, vy, vz
        DO j=1, 3
          READ(uout) tmp_dbl
          DO k=1, npart
            !dblvar(part_ind(i2)+k,j+3) = tmp_dbl(k)
            vp(part_ind(i2)+k,j) = tmp_dbl(k)
          ENDDO
        ENDDO

        !mass
        READ(uout) tmp_dbl
        DO k=1, npart
          !dblvar(part_ind(i2)+k,7) = tmp_dbl(k)
          mp(part_ind(i2)+k) = tmp_dbl(k)
        ENDDO

        !ID
        IF(larr(20) .GE. 10) THEN !long long ID
          READ(uout) tmp_l64
          DO k=1, npart
            idvar(part_ind(i2)+k) = tmp_l64(k)
          ENDDO
        ELSE !long ID
          READ(uout) tmp_int
          DO k=1, npart
            idvar(part_ind(i2)+k) = tmp_int(k)
          ENDDO
        ENDIF

        !Level
        READ(uout) tmp_int

        !FAM / TAG
        IF(larr(19) .GE. 10) THEN !READ THEM DIRECTLY
          READ(uout) tmp_byte !FAM
          DO k=1, npart
            !lonvar(part_ind(i2)+k,1) = tmp_byte(k)
            fam(part_ind(i2)+k) = tmp_byte(k)
            IF(fam(part_ind(i2)+k) .GT. 100) &
                fam(part_ind(i2)+k) = fam(part_ind(i2)+k) - 255
            !IF(lonvar(part_ind(i2)+k,1) .GT. 100) &
            !    lonvar(part_ind(i2)+k,1) = &
            !            lonvar(part_ind(i2)+k,1) - 255
          ENDDO

          READ(uout) tmp_byte !TAG
          DO k=1, npart
            !lonvar(part_ind(i2)+k,2) = tmp_byte(k)
            !IF(lonvar(part_ind(i2)+k,2) .GT. 100) &
            !    lonvar(part_ind(i2)+k,2) = &
            !            lonvar(part_ind(i2)+k,2) - 255
            tag(part_ind(i2)+k) = tmp_byte(k)
            IF(tag(part_ind(i2)+k) .GT. 100) &
                tag(part_ind(i2)+k) = tag(part_ind(i2)+k) - 255
          ENDDO
        ENDIF

        !AGE
        IF(larr(17) .LT. 10) THEN
          READ(uout) tmp_dbl
          DO k=1, npart
            !dblvar(part_ind(i2)+k,8) = tmp_dbl(k)
            ap(part_ind(i2)+k) = tmp_dbl(k)
          ENDDO
        ENDIF

        !METALLICITY
        IF(larr(18) .LT. 10) THEN
          READ(uout) tmp_dbl
          DO k=1, npart
            !dblvar(part_ind(i2)+k,9) = tmp_dbl(k)
            zp(part_ind(i2)+k) = tmp_dbl(k)
          ENDDO
        ENDIF

        !DOMAIN
        IF(larr(16) .LT. 10) THEN
        DO k=1, npart
          !lonvar(part_ind(i2)+k,3) = i
          domain(part_ind(i2)+k) = i
        ENDDO
        ENDIF

        !FAM Adding
        IF(larr(19) .LT. 10) THEN

          DO k=1, npart
          !DM
          !IF(dblvar(part_ind(i2)+k,8) .EQ. 0.) THEN
          IF(ap(part_ind(i2)+k) .EQ. 0 .AND. &
                idvar(part_ind(i2)+k).GT.0)THEN
            !lonvar(part_ind(i2)+k,1) = 1
            fam(part_ind(i2)+k) = 1
          ENDIF

          IF(ap(part_ind(i2)+k) .LT. 0.) THEN
            !lonvar(part_ind(i2)+k,1) = 2
            fam(part_ind(i2)+k) = 2
          ENDIF
          !STAR
          
          !SINK (NEGATIVE ID)

          !TRACER

          ENDDO
        ENDIF

        DEALLOCATE(tmp_dbl)
        DEALLOCATE(tmp_l64)
        DEALLOCATE(tmp_int)
        DEALLOCATE(tmp_byte)
        CLOSE(uout)
      ENDDO
      !$OMP END PARALLEL DO



      END SUBROUTINE
