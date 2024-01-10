!234567
      SUBROUTINE jsrd_part_totnum(larr, darr, &
        fname, npart_tot, part_ind, domlist)

      USE omp_lib

      IMPLICIT NONE
      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) darr(20)

      CHARACTER*(larr(4)) fname
      INTEGER(KIND=4) npart_tot
      INTEGER(KIND=4) part_ind(larr(1))
      INTEGER(KIND=4) domlist(larr(1))
!!!!!
!! LOCAL VARIABLES
!!!!!
      INTEGER(KIND=4) i, j, k, ncpu
      INTEGER(KIND=4) npartp, n_thread, cpu0, cpu1, uout
      INTEGER(KIND=4) parts(larr(1))
      CHARACTER*(100) domnum, fdum
      ncpu      = larr(1)
      n_thread  = larr(3)

      CALL OMP_SET_NUM_THREADS(n_thread)

      !$OMP PARALLEL DO default(shared) &
      !$OMP & schedule(static) private(uout,fdum,domnum,npartp,i) &
      !$OMP & reduction(+:npart_tot)
      DO j=1, ncpu
        i = domlist(j)
        WRITE(domnum,'(I5.5)') i
        uout =  OMP_GET_THREAD_NUM()
        uout = uout + 10
        fdum = TRIM(fname)//TRIM(domnum)
        OPEN(NEWUNIT=uout, FILE=fdum, FORM='unformatted', STATUS='OLD')
        READ(uout); READ(uout); READ(uout) npartp
        CLOSE(uout)
        parts(j) = npartp
        npart_tot = npart_tot + npartp
      ENDDO
      !$OMP END PARALLEL DO

      part_ind(1) = 0
      IF(ncpu .GT. 1) THEN
        DO i=2, ncpu
          DO j=1, i-1
            part_ind(i) = part_ind(i) + parts(j)
          ENDDO
        ENDDO
      ENDIF

      END SUBROUTINE

