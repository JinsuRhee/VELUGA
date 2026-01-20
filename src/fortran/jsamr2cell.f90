!234567
      SUBROUTINE jsamr2cell(larr, darr, fname_a, fname_h, fname_i, &
                    mg_ind, mesh_xg, mesh_dx, mesh_hd, mesh_lv, mesh_mp, &
                    domlist, levind, chemind, dustind)

      USE omp_lib
      IMPLICIT NONE

      INTEGER(KIND=4) larr(30)
      REAL(KIND=8) darr(20)

      CHARACTER(larr(4)) fname_a
      CHARACTER(larr(5)) fname_h
      CHARACTER(larr(6)) fname_i

      INTEGER(KIND=4) mg_ind(larr(7))
      REAL(KIND=8) mesh_xg(larr(11),larr(8))
      REAL(KIND=8) mesh_dx(larr(11))
      REAL(KIND=8) mesh_mp(larr(11))
      REAL(KIND=8) mesh_hd(larr(11),larr(12))
      INTEGER(KIND=4) mesh_lv(larr(11))
      INTEGER(KIND=4) domlist(larr(1))
      INTEGER(KIND=4) levind(larr(10),3)
      INTEGER(KIND=4) chemind(larr(22)), dustind(larr(24))

!!!!! LOCAL VARIABLES

      INTEGER(KIND=4) i, j, k, l, ndom, n_thread, ncpu, ndim, levelmin, levelmax
      INTEGER(KIND=4) nx, ny, nz, nvarh, ntot, nboundary, nvarhnew
      INTEGER(KIND=4) icpu, omp_ind, levind_tmp(larr(10))

      INTEGER(KIND=4) ilev, i0, i1, ind0, ind1
      INTEGER(KIND=4) levind2(larr(1), larr(10))
      INTEGER(KIND=4) lind_1dtbl(larr(10)*larr(1),2), lind_maxind
      INTEGER(KIND=4) ncell_onlevel(larr(10)), ind_onlevel(larr(10),2)
      INTEGER(KIND=4) omp_levind(larr(10))

      REAL(KIND=8) mesh_xg2(larr(11),larr(8))
      REAL(KIND=8) mesh_dx2(larr(11))
      REAL(KIND=8) mesh_hd2(larr(11),larr(12))
      INTEGER(KIND=4) mesh_lv2(larr(11))
      REAL(KIND=8) tokpc, tokms, toKmu, tocc, tomsun
      INTEGER(KIND=4) skipchem, skipdust, nvarnew, ivar, ivar0
      INTEGER(KIND=4) nchem, ndust
      LOGICAL ok_ind

      ndom        = larr(1)
      n_thread    = larr(3)
      ncpu        = larr(7)
      ndim        = larr(8)
      levelmin    = larr(9)
      levelmax    = larr(10)
      ntot        = larr(11)
      nvarh       = larr(12)
      nx          = larr(13)
      ny          = larr(14)
      nz          = larr(15)
      nboundary   = larr(16)
      nvarhnew    = larr(17)
      skipchem    = larr(21)
      nchem       = larr(22)
      skipdust    = larr(23)
      ndust       = larr(24)

      tokpc       = darr(1)
      tokms       = darr(2)
      tomsun      = darr(3)
      tocc        = darr(4)
      toKmu       = darr(5)


      CALL OMP_SET_NUM_THREADS(n_thread)
 

      !$OMP PARALLEL DO &
      !$OMP & shared(larr, fname_a, fname_h, ndom) &
      !$OMP & shared(mesh_xg, mesh_dx, mesh_lv, mesh_hd, levind2) &
      !$OMP & private(icpu, omp_ind, levind_tmp)
      DO i=1, ndom
        icpu     = domlist(i)
        omp_ind  = mg_ind(domlist(i))

        levind_tmp  = 0

        CALL jsamr2cell_read(larr, icpu, omp_ind, fname_a, fname_h, &
                mesh_xg2, mesh_dx2, mesh_lv2, mesh_hd2, levind_tmp, &
                chemind, dustind)

        levind2(i,:)   = levind_tmp
      ENDDO
      !$OMP END PARALLEL DO


      !!-----
      !! SORTING BY LEVEL
      !!-----

      !!----- First Count cell on each level

      lind_1dtbl = 0
      ind0  = 1
      DO i=1, ndom
        DO ilev=levelmin, levelmax
          IF(levind2(i,ilev) .EQ. 0) CYCLE

          lind_1dtbl(ind0,1) = levind2(i,ilev)
          lind_1dtbl(ind0,2) = ilev
          ind0 = ind0 + 1
        ENDDO
      ENDDO
      lind_maxind = ind0-1
      
      ncell_onlevel = 0
      DO i=1, lind_maxind
        
        IF(i.EQ.1) THEN
          i0 = 0
        ELSE
          i0 = lind_1dtbl(i-1,1)
        ENDIF

        ncell_onlevel(lind_1dtbl(i,2)) = ncell_onlevel(lind_1dtbl(i,2)) + lind_1dtbl(i,1)-i0
      ENDDO

      ind0 = 0
      DO i=1, levelmax
        ind_onlevel(i,1) = ind0+1
        ind0 = ind0 + ncell_onlevel(i)
        ind_onlevel(i,2) = ind0

        levind(i,1) = ind_onlevel(i,1)-1
        levind(i,2) = ind_onlevel(i,2)-1
        levind(i,3) = ncell_onlevel(i)
      ENDDO

      !!----- ORDERING
      

      !shared ind_onlevel, levind2, mesh_xg, mesh_xg2, mesh...
      !private ind0, ind1, i, i0, i1
      !$OMP PARALLEL DO &
      !$OMP & shared(ind_onlevel, levind2, ndom, levelmin, levelmax, nvarh) &
      !$OMP & shared(mesh_xg, mesh_dx, mesh_lv, mesh_hd, mesh_mp) &
      !$OMP & shared(mesh_xg2, mesh_dx2, mesh_lv2, mesh_hd2) &
      !$OMP & shared(nchem, ndust, nvarhnew) &
      !$OMP & private(ind0, ind1, i, i0, i1, j, l, ivar0, ivar, ok_ind)
      DO ilev=levelmin, levelmax
        ind0  = ind_onlevel(ilev,1) !! starting index


        DO i=1, ndom
          IF(levind2(i,ilev) .EQ. 0) CYCLE

          !! INDICES ON CELL ARRAY
          i0 = levind2(i,ilev-1)
          i1 = levind2(i,ilev)
          IF(i0 .EQ. 0 .AND. i .NE. 1) THEN
            DO k=i-1,1,-1
              DO j=levelmax, levelmin, -1
                IF(levind2(k,j) .NE. 0) THEN
                  i0 = levind2(k,j)
                  EXIT
                ENDIF
              ENDDO
              IF(i0 .NE. 0) EXIT
            ENDDO
          ENDIF
          i0 = i0 + 1

          ind1  = i1-i0 + ind0
          !print *, i, ilev, i0, i1, ind0, ind1, size(mesh_lv2)

          mesh_xg(ind0:ind1,:)  = mesh_xg2(i0:i1,:)*tokpc
          mesh_lv(ind0:ind1)    = mesh_lv2(i0:i1)
          mesh_dx(ind0:ind1)    = mesh_dx2(i0:i1)*tokpc

          mesh_hd(ind0:ind1,2:4)= mesh_hd2(i0:i1,2:4)*tokms
          mesh_hd(ind0:ind1,1)  = mesh_hd2(i0:i1,1)*tocc
          mesh_hd(ind0:ind1,5)  = mesh_hd2(i0:i1,5)/mesh_hd2(i0:i1,1)*toKmu
          mesh_hd(ind0:ind1,6)  = mesh_hd2(i0:i1,6)

          mesh_mp(ind0:ind1)    = mesh_hd2(i0:i1,1)*(mesh_dx2(i0:i1)**3.)*tomsun

          IF(nvarh .GE. 7) THEN
            ivar0 = 7
            DO ivar=7, nvarh
              ok_ind = .true.
              IF(skipchem .EQ. 1) THEN
                DO l=1, nchem
                  IF(ivar .EQ. chemind(l)+1) THEN
                    ok_ind = .false.
                    EXIT
                  ENDIF
                ENDDO
              ENDIF

              IF(skipdust .EQ. 1) THEN
                DO l=1, ndust
                  IF(ivar .EQ. dustind(l)+1) THEN
                    ok_ind = .false.
                    EXIT
                  ENDIF
                ENDDO
              ENDIF

              IF(ok_ind) THEN
                mesh_hd(ind0:ind1,ivar0) = mesh_hd2(i0:i1,ivar0)
                ivar0 = ivar0 + 1
              ENDIF
            ENDDO

          ENDIF

          ind0 = ind1 + 1
        ENDDO
      ENDDO
      END SUBROUTINE

      !!-----
      !! READING ROUTINE
      !!-----
      SUBROUTINE jsamr2cell_read(larr, icpu, omp_ind, fname_a, fname_h, &
        mesh_xg, mesh_dx, mesh_lv, mesh_hd, levind_tmp, chemind, dustind)

      USE OMP_lib
      IMPLICIT NONE
      
      INTEGER(KIND=4) larr(30), icpu, omp_ind

      REAL(KIND=8) mesh_xg(larr(11),larr(8))
      REAL(KIND=8) mesh_dx(larr(11))
      REAL(KIND=8) mesh_hd(larr(11),larr(12))
      INTEGER(KIND=4) mesh_lv(larr(11))
      INTEGER(KIND=4) levind_tmp(larr(10))
      INTEGER(KIND=4) chemind(larr(22)), dustind(larr(24))

      

      CHARACTER(larr(4)) fname_a
      CHARACTER(larr(5)) fname_h

      !!-----
      !! LOCAL VARIABLES
      !!-----
      INTEGER(KIND=4) ncpu, nboundary, levelmax, ndim, nvarh, ivar
      INTEGER(KIND=4) i, j, k, l, m, ilevel, ind, idim
      INTEGER(KIND=4) uout, uout2
      INTEGER(KIND=4) ngrid(1:larr(7)+larr(16),1:larr(10))
      INTEGER(KIND=4) ngrid_b(1:larr(16),1:larr(10))
      INTEGER(KIND=4) nx_full, ny_full, nz_full

      INTEGER(KIND=4) ix, iy, iz, ngrida, nx, ny, nz, twotondim
      INTEGER(KIND=4) force_levcut
      INTEGER(KIND=4) skipchem, skipdust, nvarhnew
      INTEGER(KIND=4) ivar0

      REAL(KIND=8), DIMENSION(1:3) :: xbound=(/0d0,0d0,0d0/)
      REAL(KIND=8), DIMENSION(1:8,1:3) :: xc
      REAL(KIND=8) dx, dx2
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::rho
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::x, xg
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE ::var
      INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: son
      LOGICAL ok_cell, ok_ind
      

      INTEGER(KIND=4) nchem, ndust

      CHARACTER(100) domnum, fdum_a, fdum_h, ordering

      !!----- Initial set
      !icpu      = larr(17)
      ncpu      = larr(7)
      ndim      = larr(8)
      nvarh     = larr(12)
      !omp_ind   = larr(18)
      nboundary = larr(16)
      levelmax  = larr(10)
      nx        = larr(13)
      ny        = larr(14)
      nz        = larr(15)
      nvarhnew  = larr(17)

      force_levcut = larr(20)
      
      skipchem    = larr(21)
      nchem       = larr(22)
      skipdust    = larr(23)
      ndust       = larr(24)

      ngrid   = 0
      ordering = 'hilbert'
      twotondim = 2**ndim

      !!----- OPEN FILEs
      WRITE(domnum, '(I5.5)') icpu

      fdum_a = TRIM(fname_a)//TRIM(domnum)
      fdum_h = TRIM(fname_h)//TRIM(domnum)


      uout = OMP_GET_THREAD_NUM() + 10
      uout2 = OMP_GET_THREAD_NUM() + 10 + ncpu

      OPEN(UNIT=uout,  FILE=fdum_a, FORM='unformatted', STATUS='old')
      OPEN(UNIT=uout2, FILE=fdum_h, FORM='unformatted', STATUS='old')
      
      !!----- READ first lines of AMR
      DO i=1, 21
        READ(uout)
      ENDDO

      !!----- Read ngrids
      READ(uout) ngrid(1:ncpu,1:levelmax)

      READ(uout)
      IF(nboundary>0)THEN
        READ(uout); READ(uout)
        READ(uout) ngrid_b
        ngrid(ncpu+1:ncpu+nboundary,1:levelmax)=ngrid_b
      ENDIF

      READ(uout); READ(uout)


      !!----- Do not take care of bisection ordering here
      IF(TRIM(ordering) .EQ. 'bisection') THEN
        READ(uout); READ(uout); READ(uout); READ(uout); READ(uout)
      ELSE
        READ(uout)
      ENDIF
      READ(uout); READ(uout); READ(uout)

      
      !!----- READ Hydro before going to level loop
      READ(uout2); READ(uout2); READ(uout2)
      READ(uout2); READ(uout2); READ(uout2)

      !!----- Loop Over levels
      DO ilevel=1, levelmax
        dx = 0.5 ** ilevel
        dx2 = 0.5*dx

        nx_full = 2**ilevel
        ny_full = 2**ilevel
        nz_full = 2**ilevel


        DO ind=1, twotondim
          iz = (ind-1)/4
          iy = (ind-1-4*iz)/2
          ix = (ind-1-2*iy-4*iz)
          xc(ind,1) = (dble(ix)-0.5D0)*dx 
          xc(ind,2) = (dble(iy)-0.5D0)*dx
          xc(ind,3) = (dble(iz)-0.5D0)*dx
        ENDDO

        ngrida = ngrid(icpu,ilevel)

        !!----- ALLOCATE
        IF(ngrida>0) THEN
          ALLOCATE(xg (1:ngrida,1:ndim))
          ALLOCATE(son(1:ngrida,1:twotondim))
          ALLOCATE(var(1:ngrida,1:twotondim,1:nvarh))
          ALLOCATE(x  (1:ngrida,1:ndim))
          ALLOCATE(rho(1:ngrida))
        ENDIF

        !!----- Loop Over domains
        DO j=1,nboundary+ncpu
          IF(ngrid(j,ilevel)>0) THEN
            READ(uout); READ(uout); READ(uout)
            DO idim=1, ndim !READ GRID CENTER
              IF(j.EQ.icpu)THEN
                READ(uout) xg(:,idim)
              ELSE
                READ(uout)
              ENDIF
            ENDDO

            READ(uout) ! SKIP FATHER INDEX

            DO ind=1, 2*ndim ! SKIP NBOR INDEX
              READ(uout)
            ENDDO

            DO ind=1, twotondim ! READ SON INDEX
              IF(j.EQ.icpu)THEN
                READ(uout) son(:,ind)
              ELSE
                READ(uout)
              ENDIF
            ENDDO

            DO ind=1, twotondim ! SKIP CPU MAP
              READ(uout)
            ENDDO

            DO ind=1, twotondim !SKIP REFINEMENT MAP
              READ(uout)
            ENDDO
          ENDIF

          !!READ HYDRO
          READ(uout2); READ(uout2)
            
          IF(ngrid(j,ilevel)>0)THEN
            DO ind=1, twotondim
              DO ivar=1, nvarh
                IF(j.EQ.icpu)THEN
                  READ(uout2) var(:,ind,ivar)
                ELSE
                  READ(uout2)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO

        !!----- Loop Over cells
        IF(ngrida>0) THEN
          DO k=1, ngrida !MERGE DATA
            DO ind=1, twotondim
              x(k,1) = (xg(k,1)+xc(ind,1)-xbound(1))
              x(k,2) = (xg(k,2)+xc(ind,2)-xbound(2))
              x(k,3) = (xg(k,3)+xc(ind,3)-xbound(3))

              IF(force_levcut .LT.0) THEN
                ok_cell = .NOT. (son(k,ind)>0 .AND. ilevel<levelmax)
              ELSE
                ok_cell = (ilevel<=force_levcut) .AND. (son(k,ind)==0 .OR. ilevel==force_levcut)
              ENDIF
              IF(ok_cell) THEN
                omp_ind = omp_ind + 1

                mesh_xg(omp_ind,1) = x(k,1)
                mesh_xg(omp_ind,2) = x(k,2)
                mesh_xg(omp_ind,3) = x(k,3)
                mesh_dx(omp_ind)   = dx
                mesh_lv(omp_ind)   = ilevel

                levind_tmp(ilevel)  = omp_ind
                ivar0 = 1
                DO ivar=1, nvarh
                  ok_ind = .true.
                  IF(skipchem .EQ. 1) THEN
                    DO l=1, nchem
                      IF(ivar .EQ. chemind(l)+1) THEN
                        ok_ind = .false.
                        EXIT
                      ENDIF
                    ENDDO
                  ENDIF

                  IF(skipdust .EQ. 1) THEN
                    DO l=1, ndust
                      IF(ivar .EQ. dustind(l)+1) THEN
                        ok_ind = .false.
                        EXIT
                      ENDIF
                    ENDDO
                  ENDIF

                  IF(ok_ind) THEN
                    mesh_hd(omp_ind,ivar0) = var(k, ind, ivar)
                    ivar0 = ivar0 + 1
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        
          DEALLOCATE(xg, son, var, rho, x)
        ENDIF
      ENDDO

      CLOSE(uout)
      CLOSE(uout2)

      END SUBROUTINE
