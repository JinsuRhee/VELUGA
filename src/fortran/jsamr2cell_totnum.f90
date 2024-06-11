!234567
      SUBROUTINE jsamr2cell_totnum(larr, darr, fname_a, fname_h, &
                      ntot, nvarh, mg_ind, domlist)

      USE omp_lib
      IMPLICIT NONE

      INTEGER(KIND=4) larr(20)
      REAL(KIND=8) darr(20)

      CHARACTER(larr(4)) fname_a
      CHARACTER(larr(5)) fname_h

      INTEGER(KIND=4) ntot
      INTEGER(KIND=4) nvarh
      INTEGER(KIND=4) mg_ind(larr(6)), domlist(larr(1))

!!!!! LOCAL VARIABLES
      INTEGER(KIND=4) i, j, icpu, k, idim, ind, nx, ny, nz
      INTEGER(KIND=4) n_thread, cpu0, cpu1, uout, uout2

      !!INFO RELATED
      INTEGER(KIND=4) ncpu, ndim, levelmin, levelmax, ndom, nboundary
      !!AMR RELATED
      INTEGER(KIND=4) mg_num(larr(6),larr(9))
      INTEGER(KIND=4) mg_indtmp(larr(6))
      INTEGER(KIND=4) ind_tmp, ngrida, twotondim, ntemp

      INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridfile
      INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ngridbound
      INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: son

      CHARACTER(100) domnum, fdum_a, fdum_h

      ndom = larr(1)
      n_thread = larr(3)
      ncpu = larr(6)
      ndim = larr(7)
      levelmin = larr(8)
      levelmax = larr(9)
      twotondim = 2**ndim
      CALL OMP_SET_NUM_THREADS(n_thread)

      mg_indtmp = 0
      !!!!!
      !! READ HYDRO/AMR HEADERS
      !!!!!
      WRITE(domnum, '(I5.5)') domlist(1)
      fdum_h = TRIM(fname_h)//TRIM(domnum)
      OPEN(UNIT=10, FILE=fdum_h, FORM='unformatted', STATUS='old')
      READ(10); READ(10) nvarh; CLOSE(10)

      fdum_a  = TRIM(fname_a)//TRIM(domnum)
      OPEN(UNIT=11, FILE=fdum_a, FORM='unformatted', STATUS='old')
      READ(11); READ(11); READ(11) nx, ny, nz
      READ(11); READ(11); READ(11) nboundary
      READ(11); READ(11)
      CLOSE(11)


      larr(11)  = nx
      larr(12)  = ny
      larr(13)  = nz
      larr(14)  = nboundary


      IF(nboundary>0)ALLOCATE(ngridbound(1:nboundary,1:levelmax))
      ALLOCATE(ngridfile(1:ncpu+nboundary,1:levelmax))

      ntot = 0
      !!!!!
      !! COUNT TOTAL NUM
      !!!!!
      !$OMP PARALLEL DO default(shared) &
      !$OMP & schedule(static) private(domnum, fdum_a) &
      !$OMP & private(uout, mg_num, ind_tmp, i, j, k, icpu, idim) &
      !$OMP & private(ngridfile, ngridbound, ngrida) &
      !$OMP & reduction(+: ntot)
      DO i=1, ndom
        icpu = domlist(i)

        WRITE(domnum, '(I5.5)') icpu
        fdum_a = TRIM(fname_a)//TRIM(domnum)
        uout = OMP_GET_THREAD_NUM() + 10
        OPEN(UNIT=uout, FILE=fdum_a, FORM='unformatted', STATUS='old')

        DO j=1, 21
          READ(uout)
        ENDDO

        READ(uout) ngridfile(1:ncpu,1:levelmax)

        READ(uout)
        IF(nboundary>0)THEN
          READ(uout); READ(uout)
          READ(uout) ngridbound
          ngridfile(ncpu+1:ncpu+nboundary,1:levelmax)=ngridbound
        ENDIF
        READ(uout); READ(uout)

        ! Do not take care of bisection ordering here
        READ(uout); READ(uout); READ(uout); READ(uout)

        

        ! Loop over levels
        DO j=1, levelmax
          ngrida  = ngridfile(icpu,j)
          IF(ngrida>0) THEN
            ALLOCATE(son(1:ngrida,1:twotondim))
          ENDIF


          ! Loop over domains
          DO k=1,nboundary+ncpu
            IF(ngridfile(k,j)>0)THEN
              READ(uout); READ(uout); READ(uout)
              DO idim=1, ndim
                ! SKIP READ GRID CENTER
                READ(uout)
              ENDDO

              READ(uout) ! SKIP FATHER INDEX

              DO ind=1, 2*ndim ! SKIP NBOR INDEX
                READ(uout)
              ENDDO

              DO ind=1, twotondim ! READ SON INDEX
                IF(k.EQ.icpu)THEN
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
          ENDDO

          IF(ngrida>0) THEN
            DO k=1, ngrida !MERGE DATA

              DO ind=1, twotondim
                IF(son(k,ind)==0) THEN
                  ntot = ntot + 1
                  mg_indtmp(icpu) = mg_indtmp(icpu) + 1
                ENDIF
              ENDDO
            ENDDO
          ENDIF

          IF(ngrida>0) THEN
            DEALLOCATE(son)
          ENDIF

        ENDDO
        CLOSE(uout)
      ENDDO
      !$OMP END PARALLEL DO


      ind_tmp   = 0
      DO i=1, ndom
        icpu = domlist(i)
  
        mg_ind(icpu) = ind_tmp
        ind_tmp = ind_tmp + mg_indtmp(icpu)
      ENDDO

      IF(nboundary>0)DEALLOCATE(ngridbound)
      DEALLOCATE(ngridfile)
      END SUBROUTINE
