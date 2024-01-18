MODULE read_ramses_py
    !REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: p_sgl
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: r_dbl
    !INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: r_int
    INTEGER(KIND=8), DIMENSION(:,:), ALLOCATABLE :: r_lnt

    TYPE ramses_type
        INTEGER(KIND=4) idtype      !! 1 (long) -1 (long64)
        INTEGER(KIND=4) famtype     !! 1 (new) -1 (old)
        INTEGER(KIND=4) parttype    !! -1 (all) 1 (DM) 2 (STAR) 3 (sink) 

        INTEGER(KIND=4) skip_domain !! Skip domain (if a value is positive)
        INTEGER(KIND=4) skip_time   !! Skip time (if raw data do not have)
        INTEGER(KIND=4) skip_metal  !! Skip metal

        CHARACTER(LEN=1000) :: dir_raw

        REAL(KIND=8) dmp_mass
        INTEGER(KIND=4) n_snap
    END TYPE ramses_type
CONTAINS

!!-----
!! READ RAW PARTICLE
!!-----
SUBROUTINE read_part(ftype, dom_list, n_thread)
    USE omp_lib
    IMPLICIT NONE

    INTEGER(KIND=4), DIMENSION(:), INTENT(IN) :: dom_list
    INTEGER(KIND=4), INTENT(IN) :: n_thread
    TYPE(ramses_type), INTENT(IN) :: ftype

    !!-----
    !! LOCAL VARIABLES
    !!-----
    INTEGER(KIND=4) i, j, k
    INTEGER(KIND=4) n_raw, n_mpi, nn, nbody
    INTEGER(KIND=4) uout
    CHARACTER*(100) fname, domnum
    CHARACTER(LEN=5) snap
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: dom_ind, nbodyarr

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: dum_dbl
    INTEGER(KIND=8), ALLOCATABLE, DIMENSION(:) :: dum_lint
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: dum_int
    INTEGER(KIND=1), ALLOCATABLE, DIMENSION(:) :: dum_byte

    CALL OMP_SET_NUM_THREADS(n_thread)

    n_mpi = SIZE(dom_list)

    !!-----
    !! GET N_RAW
    !!-----
    ALLOCATE(dom_ind(1:n_mpi))
    ALLOCATE(nbodyarr(1:n_mpi))
    dom_ind = 0
    nbodyarr = 0

    n_raw = 0
    DO i=1, n_mpi
        WRITE(snap, '(I5.5)') ftype%n_snap
        WRITE(domnum, '(I5.5)') dom_list(i)
        fname = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/part_'// &
            TRIM(snap)//'.out'//TRIM(domnum)

        OPEN(newunit=uout, file=fname, form='unformatted', status='old')
        READ(uout); READ(uout); READ(uout) nbody
        CLOSE(uout)

        nbodyarr(i) = nbody
        n_raw = n_raw + nbody

        
    ENDDO

    IF(n_mpi .GT. 1) THEN
        DO i=2, n_mpi
            DO j=1, i-1
                dom_ind(i) = dom_ind(i) + nbodyarr(j)
            ENDDO
        ENDDO
    ENDIF
    

    DEALLOCATE(nbodyarr)
    !!-----
    !! READ RAW DATA
    !!-----
    CALL allocate_dbl(n_raw, 9)
    CALL allocate_lint(n_raw, 3)

    nn = 0

    !$OMP PARALLEL DO default(shared) &
    !$OMP & schedule(static) &
    !$OMP & private(snap, domnum, fname, uout) &
    !$OMP & private(dum_dbl, dum_int, dum_lint) &
    !$OMP & private(dum_byte, i, j, nbody)
    DO k=1, n_mpi
        WRITE(snap, '(I5.5)') ftype%n_snap
        WRITE(domnum, '(I5.5)') dom_list(k)
        fname = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/part_'// &
            TRIM(snap)//'.out'//TRIM(domnum)


        uout = OMP_GET_THREAD_NUM() + 10

        OPEN(unit=uout, file=fname, form='unformatted', status='old')
        READ(uout); READ(uout); READ(uout) nbody; READ(uout);
        READ(uout); READ(uout); READ(uout); READ(uout)

        ALLOCATE(dum_dbl(1:nbody))
        ALLOCATE(dum_lint(1:nbody))
        ALLOCATE(dum_int(1:nbody))
        ALLOCATE(dum_byte(1:nbody))


        !! Position
        DO j=1, 3
            READ(uout) dum_dbl
            DO i=1, nbody
                r_dbl(dom_ind(k)+i,j) = dum_dbl(i)
            ENDDO
        ENDDO

        !! Velocity
        DO j=1, 3
            READ(uout) dum_dbl
            DO i=1, nbody
                r_dbl(dom_ind(k)+i,j+3) = dum_dbl(i)
            ENDDO
        ENDDO

        !! Mass
        READ(uout) dum_dbl
        DO i=1, nbody
            r_dbl(dom_ind(k)+i,7) = dum_dbl(i)
        ENDDO

        
        !! ID
        IF(ftype%idtype .GT. 0) THEN
            READ(uout) dum_lint
            DO i=1, nbody
                r_lnt(dom_ind(k)+i,1) = dum_lint(i)
            ENDDO
        ELSE IF(ftype%idtype .LT. 0) THEN
            READ(uout) dum_int
            DO i=1, nbody
                r_lnt(dom_ind(k)+i,1) = dum_int(i)
            ENDDO
        ENDIF

        
        !! LEVEL
        READ(uout) dum_int

        !! (FAM ver)
        IF(ftype%famtype .EQ. 1) THEN

            !! TAG
            READ(uout) dum_byte

            !! FAMILY
            READ(uout) dum_byte
            DO i=1, nbody
                r_lnt(dom_ind(k)+i,2) = dum_byte(i)
                IF(r_lnt(dom_ind(k)+i,2) .GT. 100) &
                    r_lnt(dom_ind(k)+i,2) = r_lnt(dom_ind(k)+i,2) - 255
            ENDDO
        ENDIF

        IF(ftype%skip_time .LT. 0) THEN
            !! AGE
            READ(uout) dum_dbl
            DO i=1, nbody
                r_dbl(dom_ind(k)+i,8) = dum_dbl(i)
            ENDDO
        ENDIF

        IF(ftype%skip_metal .LT. 0) THEN
            !! AGE
            READ(uout) dum_dbl
            DO i=1, nbody
                r_dbl(dom_ind(k)+i,9) = dum_dbl(i)
            ENDDO
        ENDIF
        !! OTHER?

        DEALLOCATE(dum_dbl, dum_lint)
        DEALLOCATE(dum_int, dum_byte)
        CLOSE(uout)
    ENDDO
    !$OMP END PARALLEL DO

    DEALLOCATE(dom_ind)

    RETURN
END SUBROUTINE


!!-----
!! ROUTINES FOR MEMORY
!!-----
SUBROUTINE allocate_dbl(npart, ndim)
    IMPLICIT NONE
    INTEGER(KIND=4) npart, ndim
    IF(.NOT. ALLOCATED(r_dbl)) ALLOCATE(r_dbl(1:npart,1:ndim))
END SUBROUTINE

SUBROUTINE deallocate_dbl()
    IMPLICIT NONE
    IF(ALLOCATED(r_dbl)) DEALLOCATE(r_dbl)
END SUBROUTINE

!SUBROUTINE allocate_int(npart, ndim)
!    IMPLICIT NONE
!    INTEGER(KIND=4) npart, ndim
!
!    IF(.NOT. ALLOCATED(r_int)) ALLOCATE(r_int(1:npart,1:ndim))
!END SUBROUTINE

!SUBROUTINE deallocate_int()
!    IMPLICIT NONE
!
!    IF(ALLOCATED(r_dbl)) DEALLOCATE(r_int)
!END SUBROUTINE

SUBROUTINE allocate_lint(npart, ndim)
    IMPLICIT NONE
    INTEGER(KIND=4) npart, ndim

    IF(.NOT. ALLOCATED(r_lnt)) ALLOCATE(r_lnt(1:npart,1:ndim))
END SUBROUTINE

SUBROUTINE deallocate_lint()
    IMPLICIT NONE

    IF(ALLOCATED(r_lnt)) DEALLOCATE(r_lnt)
END SUBROUTINE

END MODULE