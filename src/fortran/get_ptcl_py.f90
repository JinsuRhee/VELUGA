!MODULE
MODULE get_ptcl_py
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_dbl
    INTEGER(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_lnt

    INTEGER(KIND=4) n_thread, r_idtype, r_famtype, r_ptype
    INTEGER(KIND=4) r_skip_domain, r_skip_time, r_skip_metal

    REAL(KIND=8) dmp_mass
    CHARACTER(LEN=1000) dir_raw

CONTAINS

    SUBROUTINE get_ptcl(n_snap, id, dom_list)

    
    USE read_ramses_py
    USE omp_lib

    IMPLICIT NONE

    INTEGER(KIND=4), INTENT(IN) :: n_snap
    INTEGER(KIND=8), DIMENSION(:), INTENT(IN) :: id
    INTEGER(KIND=4), DIMENSION(:), INTENT(IN) :: dom_list

    CHARACTER(LEN=100) aaa

    !!-----
    !! LOCAL VARIABLES
    !!-----
    TYPE(ramses_type) ftype
    INTEGER(KIND=4) npart, n_raw
    INTEGER(KIND=8) dumint
    INTEGER(KIND=4) i, j, k, l, ind0, ind1, mi

    !!----- INITIAL SETTINGS
    ftype%idtype = r_idtype
    ftype%famtype = r_famtype
    ftype%parttype = r_ptype
    ftype%skip_domain = r_skip_domain
    ftype%skip_time = r_skip_time
    ftype%skip_metal = r_skip_metal
    ftype%dir_raw = dir_raw
    ftype%dmp_mass = dmp_mass
    ftype%n_snap = n_snap

    npart = SIZE(id)




    !!-----
    !! READ RAW PART
    !!-----
    CALL read_part(ftype, dom_list, n_thread)
    n_raw   = SIZE(r_lnt(:,1))

    !!-----
    !! ALLOCATE
    !!-----
    CALL get_ptcl_allocate(npart)


    !!-----
    !! MATCHING
    !!-----

    !!----- FIRST SORT PARTICLE BY ID
    CALL SORT_PTCL(r_dbl, r_lnt)

    !!----- MATCHING BY BINARY SEARCH
    
    CALL OMP_SET_NUM_THREADS(n_thread)
    !$OMP PARALLEL DO default(shared) &
    !$OMP & private(ind0, ind1, l, mi)
    DO i=1, npart
        IF(id(i) .LT. r_lnt(1,1) .OR. id(i) .GT. r_lnt(n_raw,1)) CYCLE
        
        ind0 = 1
        ind1 = n_raw

        DO WHILE(ind1 - ind0 .GT. 10)
            l = INT((ind1+ind0)/2)
            IF(id(i) .GT. r_lnt(l,1)) ind0 = l
            IF(id(i) .LT. r_lnt(l,1)) ind1 = l
            IF(id(i) .EQ. r_lnt(l,1)) THEN
                ind0 = l; ind1 = l; mi = l
            ENDIF
        ENDDO

        IF(ind0 .NE. ind1) THEN
            mi = -1
            DO l=ind0, ind1
                IF(id(i) .EQ. r_lnt(l,1)) THEN
                    mi = l
                    EXIT
                ENDIF
            ENDDO
        ELSE
            mi = ind0
        ENDIF

        IF(mi .GE. 1) THEN
            p_dbl(i,:) = r_dbl(mi,:)
            p_lnt(i,:) = r_lnt(mi,:)
        ENDIF
    ENDDO
    !$OMP END PARALLEL DO

    !!-----
    !! MEMORY FREE
    !!-----
    CALL deallocate_dbl()
    CALL deallocate_lint()

    END SUBROUTINE

!!-----
!! SUBROUTINES
!!-----
    SUBROUTINE SORT_PTCL(p_dbl, p_lnt)

    
    USE omp_lib
    
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:) :: p_dbl
    INTEGER(KIND=8), DIMENSION(:,:) :: p_lnt


    INTEGER(KIND=4) nptcl, i, j, i0, i1
    INTEGER(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_ind


    nptcl = SIZE(p_lnt(:,1))
    CALL OMP_SET_NUM_THREADS(n_thread)

    ALLOCATE(p_ind(1:nptcl, 1:2))
    !$OMP PARALLEL DO default(shared)
    DO i=1, nptcl
        p_ind(i,1) = p_lnt(i,1) !! input ID
        p_ind(i,2) = i          !! set index
    ENDDO
    !$OMP END PARALLEL DO

    i0 = 1
    i1 = nptcl

    CALL quicksort(p_ind, nptcl, i0, i1)

    DO j=1, SIZE(p_lnt(1,:))
        p_lnt(:,j) = p_lnt(p_ind(:,2),j)
    ENDDO


    DO j=1, SIZE(p_dbl(1,:))
        p_dbl(:,j) = p_dbl(p_ind(:,2),j)
    ENDDO

    RETURN
    END SUBROUTINE

!!!!!
!! QUICK SORT
!!!!!
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort(a, nn, first, last)
    IMPLICIT NONE
    INTEGER(KIND=8) a(nn,2), x, t, n
    INTEGER(KIND=4) first, last
    INTEGER(KIND=4) i, j, nn

    x = a( (first+last) / 2, 1)
    i = first
    j = last

    DO
        DO WHILE ( a(i,1) < x )
            i = i+1
        ENDDO

        DO WHILE ( x < a(j,1) )
            j = j-1
        ENDDO

        IF(i .GE. j) EXIT

        t = a(i,1); a(i,1) = a(j,1); a(j,1) = t
        n = a(i,2); a(i,2) = a(j,2); a(j,2) = n
        i = i+1
        j = j-1
    ENDDO
    IF (first < i-1) CALL quicksort(a, nn, first, i-1)
    IF (j+1 < last) CALL quicksort(a, nn, j+1, last)
    RETURN
    END SUBROUTINE quicksort


!!-----
!! ROUTINES FOR MEMORY
!!-----
    SUBROUTINE get_ptcl_allocate(npart)
    IMPLICIT NONE
    INTEGER(KIND=4) npart
    IF(.NOT. ALLOCATED(p_dbl)) ALLOCATE(p_dbl(1:npart,1:9))
    IF(.NOT. ALLOCATED(p_lnt)) ALLOCATE(p_lnt(1:npart, 1:3))
    END SUBROUTINE

    SUBROUTINE get_ptcl_deallocate()
    IMPLICIT NONE
    INTEGER(KIND=4) npart
    IF(ALLOCATED(p_dbl)) DEALLOCATE(p_dbl)
    IF(ALLOCATED(p_lnt)) DEALLOCATE(p_lnt)
    END SUBROUTINE

END MODULE