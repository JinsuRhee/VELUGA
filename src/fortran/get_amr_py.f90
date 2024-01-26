!MODULE
MODULE get_amr_py
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_dbl
    INTEGER(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_lnt
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: p_int

    INTEGER(KIND=4) n_thread, n_mpi, levmax, n_dim
    !INTEGER(KIND=4) r_skip_domain, r_skip_time, r_skip_metal

    !REAL(KIND=8) dmp_mass
    CHARACTER(LEN=1000) dir_raw

CONTAINS

    SUBROUTINE get_amr_box(n_snap, x0, y0, z0, dx, dom_list)

    
    USE read_ramses_py
    USE omp_lib

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN) :: n_snap
    INTEGER(KIND=4), DIMENSION(:), INTENT(IN) :: dom_list
    REAL(KIND=8), INTENT(IN) :: x0, y0, z0, dx
    !!----- LOCAL VARIABLES
    INTEGER(KIND=4) nn, i, last_ind
    TYPE(ramses_type) ftype
    LOGICAL, DIMENSION(:), ALLOCATABLE :: tag
    
    !!----- INITIAL SETTINGS
    ftype%n_mpi = n_mpi
    ftype%dir_raw = dir_raw
    ftype%n_snap = n_snap
    ftype%levmax = levmax
    ftype%n_dim = n_dim

    !!----- READ RAW CELLS

    CALL read_amr(ftype, dom_list, n_thread)

    !!-----
    !! Find specified cells within a box
    !!-----
    
    ALLOCATE(tag(1:n_cell))
    tag = .false.

    last_ind = n_dim + 1

    CALL OMP_SET_NUM_THREADS(n_thread)

    nn = 0
    !$OMP PARALLEL DO default(shared) &
    !$OMP & reduction(+:nn)
    DO i=1, n_cell
        IF(r_int(i,1) .LT. 0) CYCLE

        IF(r_dbl(i,1) .LT. x0-dx/2. - r_dbl(i,last_ind)) CYCLE
        IF(r_dbl(i,1) .GT. x0+dx/2. + r_dbl(i,last_ind)) CYCLE

        IF(r_dbl(i,2) .LT. y0-dx/2. - r_dbl(i,last_ind)) CYCLE
        IF(r_dbl(i,2) .GT. y0+dx/2. + r_dbl(i,last_ind)) CYCLE

        IF(r_dbl(i,3) .LT. z0-dx/2. - r_dbl(i,last_ind)) CYCLE
        IF(r_dbl(i,3) .GT. z0+dx/2. + r_dbl(i,last_ind)) CYCLE

        tag(i) = .true.
        nn = nn + 1
    ENDDO
    !$OMP END PARALLEL DO

    !!-----
    !! ALLOCATE
    !!-----
    CALL get_amr_allocate(nn, n_dim+1)

    nn = 1
    DO i=1, n_cell
        IF(tag(i)) THEN
            p_dbl(nn,:) = r_dbl(i,:)
            p_int(nn,:) = r_int(i,:)
            nn = nn + 1
        ENDIF
    ENDDO

    !!-----
    !! MEMORY FREE
    !!-----
    CALL deallocate_dbl()
    CALL deallocate_int()
    DEALLOCATE(tag)

    END SUBROUTINE

!!-----
!! ROUTINES FOR MEMORY
!!-----
    SUBROUTINE get_amr_allocate(npart, ndim)
    IMPLICIT NONE
    INTEGER(KIND=4) npart, ndim, nvarh
    IF(.NOT. ALLOCATED(p_dbl)) ALLOCATE(p_dbl(1:npart,1:ndim))
    IF(.NOT. ALLOCATED(p_int)) ALLOCATE(p_int(1:npart, 1))
    END SUBROUTINE

    SUBROUTINE get_amr_deallocate()
    IMPLICIT NONE
    INTEGER(KIND=4) npart
    IF(ALLOCATED(p_dbl)) DEALLOCATE(p_dbl)
    IF(ALLOCATED(p_int)) DEALLOCATE(p_int)
    END SUBROUTINE


END MODULE