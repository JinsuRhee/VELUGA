MODULE read_ramses_py
    !REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: p_sgl
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: r_dbl
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: r_int
    INTEGER(KIND=8), DIMENSION(:,:), ALLOCATABLE :: r_lnt
    INTEGER(KIND=4) n_cell, n_varh
    !! HYDRO RELATED
    !INTEGER(KIND=4) h_n, h_nvarh
    !INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: h_mg_ind

    TYPE ramses_type
        INTEGER(KIND=4) idtype      !! 1 (long) -1 (long64)
        INTEGER(KIND=4) famtype     !! 1 (new) -1 (old)
        INTEGER(KIND=4) parttype    !! -1 (all) 1 (DM) 2 (STAR) 3 (sink) 

        INTEGER(KIND=4) skip_domain !! Skip domain (if a value is positive)
        INTEGER(KIND=4) skip_time   !! Skip time (if raw data do not have)
        INTEGER(KIND=4) skip_metal  !! Skip metal

        CHARACTER(LEN=1000) :: dir_raw

        REAL(KIND=8) dmp_mass
        INTEGER(KIND=4) n_snap, n_mpi, levmax, n_dim, n_var
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
!! READ LEAF CELL
!!-----
SUBROUTINE read_cell(ftype, dom_list, n_thread)

    USE omp_lib

    IMPLICIT NONE
    INTEGER(KIND=4), DIMENSION(:), INTENT(IN) :: dom_list
    INTEGER(KIND=4), INTENT(IN) :: n_thread
    TYPE(ramses_type), INTENT(IN) :: ftype

    !!-----
    !! LOCAL VARIABLES
    !!-----
    INTEGER(KIND=4) i, i2, j, k, icpu, uout, ncpu, ivar
    INTEGER(KIND=4) ndom, ndim, ilevel, idim


    CHARACTER(1000) domnum, fdum_a, fdum_h, fdum_i, snap
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mg_num
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: mg_ind, mg_indtmp
    INTEGER(KIND=4) ind_tmp, merge_ind, merge_ind2

    !!----- INFO VAR
    CHARACTER(LEN=80) ordering
    CHARACTER(LEN=120) temp_label

    !!----- AMR VAR
    INTEGER(KIND=4) twotondim, uout_a, nx, ny, nz, ngridmax
    INTEGER(KIND=4) ngrid_current, nboundary

    REAL(KIND=8) boxlen
    REAL(KIND=8), DIMENSION(1:3) :: xbound=(/0d0,0d0,0d0/)

    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridfile
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridlevel
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridbound

    !!----- HYDRO VAR
    INTEGER(KIND=4) uout_h

    !!----- LOOP VAR
    REAL(KIND=8) dx, dx2
    INTEGER(KIND=4) nx_full, ny_full, nz_full, ind
    INTEGER(KIND=4) ix, iy, iz, ngrida
    REAL(KIND=8), DIMENSION(1:8,1:3) :: xc
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::rho
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::x, xg
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE ::var
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: son
    LOGICAL, DIMENSION(:), ALLOCATABLE ::ref
    LOGICAL ok_cell
    REAL(KIND=8) :: xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1

    REAL(KIND=8) scale
    INTEGER(KIND=4) nx_loc !, icoarse_max, icoarse_min
    !!-----
    !! INITIAL SETTINGS
    !!-----
    ndom = SIZE(dom_list)

    ALLOCATE(mg_num(1:ftype%n_mpi, 1:ftype%levmax))
    ALLOCATE(mg_indtmp(1:ftype%n_mpi))
    ALLOCATE(mg_ind(1:ftype%n_mpi))

    mg_num = 0
    mg_indtmp = 0
    mg_ind = 0
    ndim = ftype%n_dim
    ncpu = ftype%n_mpi

    !!-----
    !! READ # OF HYDRO VARIABLES
    !!-----
    WRITE(domnum, '(I5.5)') dom_list(1)
    WRITE(snap,'(I5.5)') ftype%n_snap


    fdum_h = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/hydro_'// &
        TRIM(snap)//'.out'//TRIM(domnum)

    OPEN(unit=10, file=fdum_h, form='unformatted', status='old')
    READ(10); READ(10) n_varh; CLOSE(10)

    !!-----
    !! COUNT TOTAL NUM OF CELLS
    !!-----

    CALL OMP_SET_NUM_THREADS(n_thread)

    !$OMP PARALLEL DO default(shared) &
    !$OMP & schedule(static) private(domnum, fdum_a) &
    !$OMP & private(uout, mg_num, ind_tmp, i, j, icpu)
    DO k=1, ndom
        icpu = dom_list(k)
        WRITE(domnum, '(I5.5)') icpu
        fdum_a = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/amr_'// &
            TRIM(snap)//'.out'//TRIM(domnum)

        uout = OMP_GET_THREAD_NUM() + 10
        OPEN(unit=uout, file=fdum_a, form='unformatted', status='old')

        DO j=1, 21
            READ(uout)
        ENDDO
        READ(uout) mg_num
        CLOSE(uout)

        ind_tmp = 0
        DO i=1, ncpu
            DO j=1, ftype%levmax
                ind_tmp = ind_tmp + mg_num(i,j) * (2**ndim)
            ENDDO
        ENDDO

        mg_indtmp(icpu) = ind_tmp
    ENDDO
    !$OMP END PARALLEL DO

    ind_tmp = 0

    DO i=1, ncpu
        mg_ind(i) = 0
        DO j=1, ndom
            IF(i.EQ.dom_list(j)) THEN
                mg_ind(i) = ind_tmp
                ind_tmp = ind_tmp + mg_indtmp(i)
            ENDIF
        ENDDO
    ENDDO
    n_cell = ind_tmp

    DEALLOCATE(mg_num, mg_indtmp)

    !!-----
    !! ALLOCATE
    !!  x, y, z, n_varh (rho, vx, vy, vz, T, ...) , dx
    !!-----
    CALL allocate_dbl(n_cell, ftype%n_dim + n_varh + 1)
    CALL allocate_int(n_cell, 1)

    r_int = -1
    r_dbl = 0.
    
    !!-----
    !! READ RAW CELL
    !!-----
    twotondim = 2**ndim

    !! READ INFO FIRST
    fdum_i = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/info_'// &
        TRIM(snap)//'.txt'

    OPEN(unit=10, file=fdum_i, form='formatted', status='old')
    DO i=1, 19
        READ(10,*)
    ENDDO
    READ(10, '(A14, A80)') temp_label, ordering
    CLOSE(10)

    !! LOOP FOR DOMAIN
    WRITE(snap,'(I5.5)') ftype%n_snap

    DO i2=1, ndom
        i = dom_list(i2)
        WRITE(domnum, '(I5.5)') i
        fdum_a = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/amr_'// &
        TRIM(snap)//'.out'//TRIM(domnum)
        fdum_h = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/hydro_'// &
        TRIM(snap)//'.out'//TRIM(domnum)
    
        uout_a = 10
        uout_h = 11

        icpu = i
        merge_ind = mg_ind(i)
        
        !!!! READ AMR HEADER
        OPEN(unit=uout_a, file=fdum_a, form='unformatted', status='old')
        READ(uout_a); READ(uout_a); READ(uout_a) nx, ny, nz
        READ(uout_a); READ(uout_a) ngridmax; READ(uout_a) nboundary
        READ(uout_a) ngrid_current; READ(uout_a) boxlen

        xbound = (/dble(nx/2), dble(ny/2), dble(nz/2)/)

        !nx_loc = (icoarse_max - icoarse_min)
        nx_loc=1 !! how can I extract icoarse_max, icoarse_min
        scale=boxlen/dble(nx_loc)

        ALLOCATE(ngridlevel(1:ncpu, 1:ftype%levmax))
        ALLOCATE(ngridfile(1:ncpu+nboundary,1:ftype%levmax))

        IF(nboundary>0) ALLOCATE(ngridbound(1:nboundary,1:ftype%levmax))

        DO k=1, 13
            READ(uout_a)
        ENDDO
        READ(uout_a) ngridlevel
        ngridfile(1:ncpu, 1:ftype%levmax) = ngridlevel

        READ(uout_a)
        IF(nboundary>0)THEN
          READ(uout_a); READ(uout_a)
          READ(uout_a) ngridbound
          ngridfile(ncpu+1:ncpu+nboundary,1:ftype%levmax)=ngridbound
        ENDIF
        READ(uout_a); READ(uout_a)

        IF(TRIM(ordering) .EQ. 'bisection') THEN
          READ(uout_a); READ(uout_a); READ(uout_a); READ(uout_a); READ(uout_a)
        ELSE
          READ(uout_a)
        ENDIF
        READ(uout_a); READ(uout_a); READ(uout_a)


        !!!! READ HYDRO
        OPEN(UNIT=uout_h, FILE=fdum_h, FORM='unformatted', STATUS='old')
        READ(uout_h); READ(uout_h); READ(uout_h)
        READ(uout_h); READ(uout_h); READ(uout_h)

        !!!! LOOP OVER CELLS
        DO ilevel=1, ftype%levmax
            dx = 0.5 ** ilevel
            dx2 = 0.5*dx
            nx_full = 2** ilevel
            ny_full = 2** ilevel
            nz_full = 2** ilevel

            DO ind=1, twotondim
                iz = (ind-1)/4
                iy = (ind-1-4*iz)/2
                ix = (ind-1-2*iy-4*iz)
                xc(ind,1) = (dble(ix)-0.5D0)*dx
                xc(ind,2) = (dble(iy)-0.5D0)*dx
                xc(ind,3) = (dble(iz)-0.5D0)*dx
            ENDDO

            ngrida = ngridfile(icpu, ilevel)

            !!!!!! ALLOCATE
            IF(ngrida>0) THEN
                ALLOCATE(xg (1:ngrida,1:ndim))
                ALLOCATE(son(1:ngrida,1:twotondim))
                ALLOCATE(var(1:ngrida,1:twotondim,1:n_varh))
                ALLOCATE(x  (1:ngrida,1:ndim))
                ALLOCATE(rho(1:ngrida))
                ALLOCATE(ref(1:ngrida))
            ENDIF

            

            !!!!!! LOOP OVER DOMAINS
            DO j=1,nboundary+ncpu
                !! READ AMR
                !IF(ngrida>0) THEN
                IF(ngridfile(j,ilevel)>0)THEN
                    READ(uout_a); READ(uout_a); READ(uout_a)
                    DO idim=1, ndim !READ GRID CENTER
                        IF(j.EQ.icpu)THEN
                            READ(uout_a) xg(:,idim)
                        ELSE
                            READ(uout_a)
                        ENDIF
                    ENDDO
    
                    READ(uout_a) ! SKIP FATHER INDEX
    
                    DO ind=1, 2*ndim ! SKIP NBOR INDEX
                        READ(uout_a)
                    ENDDO
    
                    DO ind=1, twotondim ! READ SON INDEX
                        IF(j.EQ.icpu)THEN
                            READ(uout_a) son(:,ind)
                        ELSE
                            READ(uout_a)
                        ENDIF
                    ENDDO
    
                    DO ind=1, twotondim ! SKIP CPU MAP
                        READ(uout_a)
                    ENDDO
    
                    DO ind=1, twotondim !SKIP REFINEMENT MAP
                        READ(uout_a)
                    ENDDO
                ENDIF

                !! HYDRO
                READ(uout_h); READ(uout_h);

                IF(ngridfile(j,ilevel)>0)THEN
                    DO ind=1, twotondim
                        DO ivar=1, n_varh
                            IF(j.EQ.icpu)THEN
                                READ(uout_h) var(:,ind,ivar)
                            ELSE
                                READ(uout_h)
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            ENDDO

            !!!!!! LOOP OVER CELLS
            IF(ngrida>0) THEN
                IF(ngrida>100) THEN
                    !$OMP PARALLEL DO default(shared) schedule(dynamic) &
                    !$OMP & private(merge_ind2, ok_cell, ind)
                    DO k=1, ngrida !MERGE DATA
                        merge_ind2 = merge_ind + (k-1)*twotondim
                        DO ind=1, twotondim
                            x(k,1) = (xg(k,1)+xc(ind,1)-xbound(1))
                            x(k,2) = (xg(k,2)+xc(ind,2)-xbound(2))
                            x(k,3) = (xg(k,3)+xc(ind,3)-xbound(3))

                            ref(k) = son(k,ind)>0 .AND. ilevel<ftype%levmax

                            merge_ind2 = merge_ind2 + 1
                            r_int(merge_ind2,1) = -1

                            ok_cell= .NOT. ref(k) .AND. &
                                & (x(k,1) + dx2) >= xmin .AND. &
                                & (x(k,2) + dx2) >= ymin .AND. &
                                & (x(k,3) + dx2) >= zmin .AND. &
                                & (x(k,1) - dx2) <= xmax .AND. &
                                & (x(k,2) - dx2) <= ymax .AND. &
                                & (x(k,3) - dx2) <= zmax

                            IF(ok_cell) THEN
                                r_dbl(merge_ind2,1) = x(k,1) * scale
                                r_dbl(merge_ind2,2) = x(k,2) * scale
                                r_dbl(merge_ind2,3) = x(k,3) * scale
                                r_dbl(merge_ind2,ftype%n_dim + n_varh + 1) = dx * scale
                                r_int(merge_ind2,1) = ilevel
                                DO ivar=1, n_varh
                                    r_dbl(merge_ind2,ftype%n_dim+ivar) = var(k, ind, ivar)
                                ENDDO
                            ENDIF
                        ENDDO
                    ENDDO
                    !$OMP END PARALLEL DO
                ELSE
                    DO k=1, ngrida !MERGE DATA
                        merge_ind2 = merge_ind + (k-1)*twotondim
                        DO ind=1, twotondim
                            x(k,1) = (xg(k,1)+xc(ind,1)-xbound(1))
                            x(k,2) = (xg(k,2)+xc(ind,2)-xbound(2))
                            x(k,3) = (xg(k,3)+xc(ind,3)-xbound(3))

                            ref(k) = son(k,ind)>0 .AND. ilevel<ftype%levmax

                            merge_ind2 = merge_ind2 + 1
                            r_int(merge_ind2,1) = -1

                            ok_cell= .NOT. ref(k) .AND. &
                                & (x(k,1) + dx2) >= xmin .AND. &
                                & (x(k,2) + dx2) >= ymin .AND. &
                                & (x(k,3) + dx2) >= zmin .AND. &
                                & (x(k,1) - dx2) <= xmax .AND. &
                                & (x(k,2) - dx2) <= ymax .AND. &
                                & (x(k,3) - dx2) <= zmax

                            IF(ok_cell) THEN
                                r_dbl(merge_ind2,1) = x(k,1) * scale
                                r_dbl(merge_ind2,2) = x(k,2) * scale
                                r_dbl(merge_ind2,3) = x(k,3) * scale
                                r_dbl(merge_ind2,ftype%n_dim + n_varh + 1) = dx * scale
                                r_int(merge_ind2,1) = ilevel
                                DO ivar=1, n_varh
                                    r_dbl(merge_ind2,ftype%n_dim+ivar) = var(k, ind, ivar)
                                ENDDO
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF ! ENDIF for omp selection

                merge_ind = merge_ind + ngrida*twotondim

                DEALLOCATE(xg, son, var, ref, rho, x)
            ENDIF ! ENDIF for ngrida>0

        ENDDO ! ENDDO for loop for levels
        CLOSE(uout_a)
        CLOSE(uout_h)

        DEALLOCATE(ngridfile, ngridlevel)
        IF(nboundary>0) DEALLOCATE(ngridbound)


    ENDDO !ENDDO for loop for CPUs
    DEALLOCATE(mg_ind)
END SUBROUTINE

!!-----
!! READ AMR
!!-----
SUBROUTINE read_amr(ftype, dom_list, n_thread)

    USE omp_lib

    IMPLICIT NONE
    INTEGER(KIND=4), DIMENSION(:), INTENT(IN) :: dom_list
    INTEGER(KIND=4), INTENT(IN) :: n_thread
    TYPE(ramses_type), INTENT(IN) :: ftype

    !!-----
    !! LOCAL VARIABLES
    !!-----
    INTEGER(KIND=4) i, i2, j, k, icpu, uout, ncpu, ivar
    INTEGER(KIND=4) ndom, ndim, ilevel, idim


    CHARACTER(1000) domnum, fdum_a, fdum_h, fdum_i, snap
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mg_num
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: mg_ind, mg_indtmp
    INTEGER(KIND=4) ind_tmp, merge_ind, merge_ind2

    !!----- INFO VAR
    CHARACTER(LEN=80) ordering
    CHARACTER(LEN=120) temp_label

    !!----- AMR VAR
    INTEGER(KIND=4) twotondim, uout_a, nx, ny, nz, ngridmax
    INTEGER(KIND=4) ngrid_current, nboundary

    REAL(KIND=8) boxlen
    REAL(KIND=8), DIMENSION(1:3) :: xbound=(/0d0,0d0,0d0/)

    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridfile
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridlevel
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridbound
    
    !!----- LOOP VAR
    REAL(KIND=8) dx, dx2
    INTEGER(KIND=4) nx_full, ny_full, nz_full, ind
    INTEGER(KIND=4) ix, iy, iz, ngrida
    REAL(KIND=8), DIMENSION(1:8,1:3) :: xc
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::rho
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::x, xg
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE ::var
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: son
    LOGICAL, DIMENSION(:), ALLOCATABLE ::ref
    LOGICAL ok_cell
    REAL(KIND=8) :: xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1

    REAL(KIND=8) scale
    INTEGER(KIND=4) nx_loc !, icoarse_max, icoarse_min 
    !!-----
    !! INITIAL SETTINGS
    !!-----
    ndom = SIZE(dom_list)

    ALLOCATE(mg_num(1:ftype%n_mpi, 1:ftype%levmax))
    ALLOCATE(mg_indtmp(1:ftype%n_mpi))
    ALLOCATE(mg_ind(1:ftype%n_mpi))

    mg_num = 0
    mg_indtmp = 0
    mg_ind = 0
    ndim = ftype%n_dim
    ncpu = ftype%n_mpi

    !!-----
    !! READ # OF HYDRO VARIABLES
    !!-----
    WRITE(domnum, '(I5.5)') dom_list(1)
    WRITE(snap,'(I5.5)') ftype%n_snap

    !!-----
    !! COUNT TOTAL NUM OF CELLS
    !!-----

    CALL OMP_SET_NUM_THREADS(n_thread)

    !$OMP PARALLEL DO default(shared) &
    !$OMP & schedule(static) private(domnum, fdum_a) &
    !$OMP & private(uout, mg_num, ind_tmp, i, j, icpu)
    DO k=1, ndom
        icpu = dom_list(k)
        WRITE(domnum, '(I5.5)') icpu
        fdum_a = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/amr_'// &
            TRIM(snap)//'.out'//TRIM(domnum)

        uout = OMP_GET_THREAD_NUM() + 10
        OPEN(unit=uout, file=fdum_a, form='unformatted', status='old')

        DO j=1, 21
            READ(uout)
        ENDDO
        READ(uout) mg_num
        CLOSE(uout)

        ind_tmp = 0
        DO i=1, ncpu
            DO j=1, ftype%levmax
                ind_tmp = ind_tmp + mg_num(i,j) * (2**ndim)
            ENDDO
        ENDDO

        mg_indtmp(icpu) = ind_tmp
    ENDDO
    !$OMP END PARALLEL DO

    ind_tmp = 0

    DO i=1, ncpu
        mg_ind(i) = 0
        DO j=1, ndom
            IF(i.EQ.dom_list(j)) THEN
                mg_ind(i) = ind_tmp
                ind_tmp = ind_tmp + mg_indtmp(i)
            ENDIF
        ENDDO
    ENDDO
    n_cell = ind_tmp

    DEALLOCATE(mg_num, mg_indtmp)

    !!-----
    !! ALLOCATE
    !!  x, y, z, dx
    !!-----
    CALL allocate_dbl(n_cell, ftype%n_dim + 1)
    CALL allocate_int(n_cell, 1)

    r_int = -1
    r_dbl = 0.
    
    !!-----
    !! READ RAW CELL
    !!-----
    twotondim = 2**ndim

    !! READ INFO FIRST
    fdum_i = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/info_'// &
        TRIM(snap)//'.txt'

    OPEN(unit=10, file=fdum_i, form='formatted', status='old')
    DO i=1, 19
        READ(10,*)
    ENDDO
    READ(10, '(A14, A80)') temp_label, ordering
    CLOSE(10)

    !! LOOP FOR DOMAIN
    WRITE(snap,'(I5.5)') ftype%n_snap

    DO i2=1, ndom
        i = dom_list(i2)
        WRITE(domnum, '(I5.5)') i
        fdum_a = TRIM(ftype%dir_raw)//'/output_'//TRIM(snap)//'/amr_'// &
        TRIM(snap)//'.out'//TRIM(domnum)
            
        uout_a = 10
        
        icpu = i
        merge_ind = mg_ind(i)
        
        !!!! READ AMR HEADER
        OPEN(unit=uout_a, file=fdum_a, form='unformatted', status='old')
        READ(uout_a); READ(uout_a); READ(uout_a) nx, ny, nz
        READ(uout_a); READ(uout_a) ngridmax; READ(uout_a) nboundary
        READ(uout_a) ngrid_current; READ(uout_a) boxlen

        xbound = (/dble(nx/2), dble(ny/2), dble(nz/2)/)

        
        !nx_loc = (icoarse_max - icoarse_min)
        nx_loc=1 !! how can I extract icoarse_max, icoarse_min
        scale=boxlen/dble(nx_loc)


        ALLOCATE(ngridlevel(1:ncpu, 1:ftype%levmax))
        ALLOCATE(ngridfile(1:ncpu+nboundary,1:ftype%levmax))

        IF(nboundary>0) ALLOCATE(ngridbound(1:nboundary,1:ftype%levmax))

        DO k=1, 13
            READ(uout_a)
        ENDDO
        READ(uout_a) ngridlevel
        ngridfile(1:ncpu, 1:ftype%levmax) = ngridlevel

        READ(uout_a)
        IF(nboundary>0)THEN
          READ(uout_a); READ(uout_a)
          READ(uout_a) ngridbound
          ngridfile(ncpu+1:ncpu+nboundary,1:ftype%levmax)=ngridbound
        ENDIF
        READ(uout_a); READ(uout_a)

        IF(TRIM(ordering) .EQ. 'bisection') THEN
          READ(uout_a); READ(uout_a); READ(uout_a); READ(uout_a); READ(uout_a)
        ELSE
          READ(uout_a)
        ENDIF
        READ(uout_a); READ(uout_a); READ(uout_a)

        !!!! LOOP OVER CELLS
        DO ilevel=1, ftype%levmax
            dx = 0.5 ** ilevel
            dx2 = 0.5*dx
            nx_full = 2** ilevel
            ny_full = 2** ilevel
            nz_full = 2** ilevel

            DO ind=1, twotondim
                iz = (ind-1)/4
                iy = (ind-1-4*iz)/2
                ix = (ind-1-2*iy-4*iz)
                xc(ind,1) = (dble(ix)-0.5D0)*dx
                xc(ind,2) = (dble(iy)-0.5D0)*dx
                xc(ind,3) = (dble(iz)-0.5D0)*dx
            ENDDO

            ngrida = ngridfile(icpu, ilevel)

            !!!!!! ALLOCATE
            IF(ngrida>0) THEN
                ALLOCATE(xg (1:ngrida,1:ndim))
                ALLOCATE(son(1:ngrida,1:twotondim))
                ALLOCATE(var(1:ngrida,1:twotondim,1:n_varh))
                ALLOCATE(x  (1:ngrida,1:ndim))
                ALLOCATE(rho(1:ngrida))
                ALLOCATE(ref(1:ngrida))
            ENDIF

            

            !!!!!! LOOP OVER DOMAINS
            DO j=1,nboundary+ncpu
                !! READ AMR
                !IF(ngrida>0) THEN
                IF(ngridfile(j,ilevel)>0)THEN
                    READ(uout_a); READ(uout_a); READ(uout_a)
                    DO idim=1, ndim !READ GRID CENTER
                        IF(j.EQ.icpu)THEN
                            READ(uout_a) xg(:,idim)
                        ELSE
                            READ(uout_a)
                        ENDIF
                    ENDDO
    
                    READ(uout_a) ! SKIP FATHER INDEX
    
                    DO ind=1, 2*ndim ! SKIP NBOR INDEX
                        READ(uout_a)
                    ENDDO
    
                    DO ind=1, twotondim ! READ SON INDEX
                        IF(j.EQ.icpu)THEN
                            READ(uout_a) son(:,ind)
                        ELSE
                            READ(uout_a)
                        ENDIF
                    ENDDO
    
                    DO ind=1, twotondim ! SKIP CPU MAP
                        READ(uout_a)
                    ENDDO
    
                    DO ind=1, twotondim !SKIP REFINEMENT MAP
                        READ(uout_a)
                    ENDDO
                ENDIF
            ENDDO

            !!!!!! LOOP OVER CELLS
            IF(ngrida>0) THEN
                IF(ngrida>100) THEN
                    !$OMP PARALLEL DO default(shared) schedule(dynamic) &
                    !$OMP & private(merge_ind2, ok_cell, ind)
                    DO k=1, ngrida !MERGE DATA
                        merge_ind2 = merge_ind + (k-1)*twotondim
                        DO ind=1, twotondim
                            x(k,1) = (xg(k,1)+xc(ind,1)-xbound(1))
                            x(k,2) = (xg(k,2)+xc(ind,2)-xbound(2))
                            x(k,3) = (xg(k,3)+xc(ind,3)-xbound(3))

                            ref(k) = son(k,ind)>0 .AND. ilevel<ftype%levmax

                            merge_ind2 = merge_ind2 + 1
                            r_int(merge_ind2,1) = -1

                            ok_cell= .NOT. ref(k) .AND. &
                                & (x(k,1) + dx2) >= xmin .AND. &
                                & (x(k,2) + dx2) >= ymin .AND. &
                                & (x(k,3) + dx2) >= zmin .AND. &
                                & (x(k,1) - dx2) <= xmax .AND. &
                                & (x(k,2) - dx2) <= ymax .AND. &
                                & (x(k,3) - dx2) <= zmax

                            IF(ok_cell) THEN
                                r_dbl(merge_ind2,1) = x(k,1)*scale
                                r_dbl(merge_ind2,2) = x(k,2)*scale
                                r_dbl(merge_ind2,3) = x(k,3)*scale
                                r_dbl(merge_ind2,ftype%n_dim + 1) = dx*scale
                                r_int(merge_ind2,1) = ilevel
                            ENDIF
                        ENDDO
                    ENDDO
                    !$OMP END PARALLEL DO
                ELSE
                    DO k=1, ngrida !MERGE DATA
                        merge_ind2 = merge_ind + (k-1)*twotondim
                        DO ind=1, twotondim
                            x(k,1) = (xg(k,1)+xc(ind,1)-xbound(1))
                            x(k,2) = (xg(k,2)+xc(ind,2)-xbound(2))
                            x(k,3) = (xg(k,3)+xc(ind,3)-xbound(3))

                            ref(k) = son(k,ind)>0 .AND. ilevel<ftype%levmax

                            merge_ind2 = merge_ind2 + 1
                            r_int(merge_ind2,1) = -1

                            ok_cell= .NOT. ref(k) .AND. &
                                & (x(k,1) + dx2) >= xmin .AND. &
                                & (x(k,2) + dx2) >= ymin .AND. &
                                & (x(k,3) + dx2) >= zmin .AND. &
                                & (x(k,1) - dx2) <= xmax .AND. &
                                & (x(k,2) - dx2) <= ymax .AND. &
                                & (x(k,3) - dx2) <= zmax

                            IF(ok_cell) THEN
                                r_dbl(merge_ind2,1) = x(k,1) * scale
                                r_dbl(merge_ind2,2) = x(k,2) * scale
                                r_dbl(merge_ind2,3) = x(k,3) * scale
                                r_dbl(merge_ind2,ftype%n_dim + 1) = dx * scale
                                r_int(merge_ind2,1) = ilevel
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF ! ENDIF for omp selection

                merge_ind = merge_ind + ngrida*twotondim

                DEALLOCATE(xg, son, var, ref, rho, x)
            ENDIF ! ENDIF for ngrida>0

        ENDDO ! ENDDO for loop for levels
        CLOSE(uout_a)

        DEALLOCATE(ngridfile, ngridlevel)
        IF(nboundary>0) DEALLOCATE(ngridbound)


    ENDDO !ENDDO for loop for CPUs
    DEALLOCATE(mg_ind)
END SUBROUTINE

!!-----
!! ROUTINES FOR MEMORY
!!-----
SUBROUTINE allocate_dbl(npart, ndim)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(in) :: npart, ndim
    IF(.NOT. ALLOCATED(r_dbl)) ALLOCATE(r_dbl(1:npart,1:ndim))
END SUBROUTINE

SUBROUTINE deallocate_dbl()
    IMPLICIT NONE
    IF(ALLOCATED(r_dbl)) DEALLOCATE(r_dbl)
END SUBROUTINE

SUBROUTINE allocate_int(npart, ndim)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(in) :: npart, ndim

    IF(.NOT. ALLOCATED(r_int)) ALLOCATE(r_int(1:npart,1:ndim))
END SUBROUTINE
SUBROUTINE deallocate_int()
    IMPLICIT NONE

    IF(ALLOCATED(r_int)) DEALLOCATE(r_int)
END SUBROUTINE

SUBROUTINE allocate_lint(npart, ndim)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: npart, ndim

    IF(.NOT. ALLOCATED(r_lnt)) ALLOCATE(r_lnt(1:npart,1:ndim))
END SUBROUTINE

SUBROUTINE deallocate_lint()
    IMPLICIT NONE

    IF(ALLOCATED(r_lnt)) DEALLOCATE(r_lnt)
END SUBROUTINE


END MODULE