!234567
      Subroutine rv_match (larr, darr, dir_raw, &
                ID, ind_b, ind_u, xp, vp, zp, ap, mp, &
                rate, dom_list)

      USE omp_lib

      Implicit none
      Integer(kind=4) larr(20)
      Real(kind=8) darr(20)

      Integer(kind=8) ID(larr(1))
      Integer(kind=4) ind_b(larr(2),2), ind_u(larr(2),2)
      Real(kind=8) xp(larr(1),3), vp(larr(1),3)
      Real(kind=8) zp(larr(1)), ap(larr(1)), mp(larr(1))
      Integer(kind=4) dom_list(larr(2),larr(5))

      REAL(kind=4) rate(larr(2))

      Character(LEN=larr(11)) dir_raw

!!!!!! Local Variables
      Integer(kind=4) i, j, k, l, gind
      Integer(kind=4) n_ptcl, n_gal, n_thread, n_split, n_mpi
      Integer(kind=4) n_snap
      Integer(kind=4) numptcl, numptcl2, numgal, numtmp
      Integer(kind=4) g_list(larr(2)), nb, np

      Real(kind=8), dimension(:,:), allocatable :: p_d, p_d2
      Integer(kind=8), dimension(:,:), allocatable :: p_i, p_i2
      INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: p_skip
      !Integer(kind=8), dimension(:), allocatable :: id_g
      !Real(kind=8), dimension(:,:), allocatable :: xp_g, vp_g
      !Real(kind=8), dimension(:), allocatable :: ap_g, zp_g, mp_g

      Integer(kind=4) dumint, ind1, ind2, n0, n1
      INTEGER(KIND=4) dom_list2(larr(2),larr(5)), matchok
      INTEGER(KIND=4) threadnum, g_num, js_order
      Real(kind=8) dmp_mass, ptype
      INTEGER(KIND=8) noptcl
      noptcl = -9223372036854775800

      n_ptcl    = larr(1)
      n_gal     = larr(2)
      n_thread  = larr(3)
      !n_split   = larr(4)
      n_mpi     = larr(5)
      n_snap    = larr(6)

      dmp_mass  = darr(1)
      ptype     = darr(2) ! 1 for star ptcls // -1 for DM ptcls

      dom_list2 = -1
      js_order = 0
      call OMP_SET_NUM_THREADS(n_thread)
      !$OMP PARALLEL DO default(shared) &
      !$OMP & private(p_skip,p_d,p_d2,p_i,p_i2,gind,k,matchok) &
      !$OMP & private(threadnum,numgal,numptcl,numptcl2,numtmp) &
      !$OMP & private(g_num,dumint,js_order,g_list,l) &
      !$OMP & schedule(dynamic,5)
      Do i=1, n_mpi
        threadnum = OMP_GET_THREAD_NUM()
        IF(threadnum .EQ. 1) THEN
                js_order = js_order+1
                !PRINT *, js_order, ' // ', n_mpi / n_thread
        ENDIF

        !!-----
        !! Get # of Galaxies in this MPI domain
        !!-----
        k=i
        CALL GAL_LIST(dom_list, n_gal, n_mpi, k, g_list, numgal)
        IF(numgal .LT. 0) CYCLE !! Skip if there is no galaxy

        !!-----
        !! Get # of Ptcls in this domain
        !!-----
        k=i
        CALL RD_PART_NBODY(dir_raw, n_snap, k, threadnum, numptcl, larr)

        !!-----
        !! Allocate Raw Ptcl arrays
        !!-----
        ALLOCATE(p_d2(1:numptcl,1:9))
        ALLOCATE(p_i2(1:numptcl,1:3))

        !!-----
        !! Read Particles
        !!-----
        k=i
        IF(larr(19) .LE. 10) CALL RD_PART(dir_raw, n_snap, k, numptcl, &
                threadnum, larr, p_d2, p_i2)
        IF(larr(19) .GT. 10) CALL RD_PART_YZiCS(dir_raw, n_snap, &
                k, numptcl, threadnum, larr, p_d2, p_i2)


        !!-----
        !! Get # of DM or Star Ptcls
        !!-----
        IF(larr(19) .LE. 10) CALL GET_PTCL_NUM(&
                p_d2(1:numptcl,7:7), p_i2(1:numptcl,2:2), &
                numptcl, numptcl2, ptype, dmp_mass)
        IF(larr(19) .GT. 10) CALL GET_PTCL_NUM_YZiCS(&
                p_d2(1:numptcl,7:7), p_d2(1:numptcl,8:8), &
                numptcl, numptcl2, ptype, dmp_mass)
        IF(numptcl2 .LE. 0) THEN
                DEALLOCATE(p_d2, p_i2)
                CYCLE !! If there is no target ptcls
        ENDIF

        !!-----
        !! Allocate Matching Ptcl Arrays
        !!-----
        ALLOCATE(p_d(1:numptcl2,1:9))
        ALLOCATE(p_i(1:numptcl2,1:3))

        !!-----
        !! Extract target (DM / Star) ptcls
        !!-----
        IF(larr(19) .LE. 10) CALL GET_PTCL(p_d, p_i, p_d2, p_i2, &
              numptcl, numptcl2, ptype, dmp_mass)
        IF(larr(19) .GT. 10) CALL GET_PTCL_YZiCS(p_d, p_i, p_d2, p_i2, &
              numptcl, numptcl2, ptype, dmp_mass)

        DEALLOCATE(p_d2)
        DEALLOCATE(p_i2)

        !!-----
        !! Sort target ptcls by their ID
        !!-----
        CALL SORT_PTCL(p_d, p_i, numptcl2)

        !!-----
        !! Match for individual galaxy
        !!-----
        DO gind=1, numgal
          g_num = g_list(gind)

          !!-----
          !! # of ptcls in this galaxy
          !!-----
          numtmp = ind_b(g_num,2) - ind_b(g_num,1) + 1
          numtmp = numtmp + ind_u(g_num,2) - ind_u(g_num,1) + 1

          !!-----
          !! Allocate ptcl memory for this galaxy
          !!-----
          ALLOCATE(p_d2(1:numtmp,1:9))
          ALLOCATE(p_i2(1:numtmp,1:1))
          ALLOCATE(p_skip(1:numtmp))

          !!-----
          !! Bring original value from the galaxy
          !!-----
          dumint = 0
          DO k=ind_b(g_num,1)+1, ind_b(g_num,2)+1
            dumint = dumint + 1
            !p_skip(dumint) = -1
            !IF(xp(k,1) .GT. -1.0D7) p_skip(dumint) = 1

            p_i2(dumint,1) = id(k)

            DO l=1, 9
              p_d2(dumint,l) = -1.0D8
            ENDDO
          ENDDO

          DO k=ind_u(g_num,1)+1, ind_u(g_num,2)+1
            dumint = dumint + 1

            p_i2(dumint,1) = id(k)

            DO l=1, 9
              p_d2(dumint,l) = -1.0D8
            ENDDO
          ENDDO

          !!-----
          !! Doing Match
          !!-----
          CALL match(p_i, p_d, p_i2, p_d2, &
                  numptcl2, numtmp, matchok, larr, darr)

          !CALL match_hash(p_i, p_d, p_i2, p_d2, &
          !        numptcl2, numtmp, matchok, noptcl)
          !!-----
          !! Write matched data into the original array
          !!-----
          dumint = 0
          DO k=ind_b(g_num,1)+1, ind_b(g_num,2)+1
            dumint = dumint + 1
            !IF(p_skip(dumint) .GE. 0) CYCLE
            IF(p_d2(dumint,1) .LT. -1.0D7) CYCLE
            xp(k,1) = p_d2(dumint,1)
            xp(k,2) = p_d2(dumint,2)
            xp(k,3) = p_d2(dumint,3)

            vp(k,1) = p_d2(dumint,4)
            vp(k,2) = p_d2(dumint,5)
            vp(k,3) = p_d2(dumint,6)

            mp(k) = p_d2(dumint,7)
            ap(k) = p_d2(dumint,8)
            zp(k) = p_d2(dumint,9)

          ENDDO

          DO k=ind_u(g_num,1)+1, ind_u(g_num,2)+1
            dumint = dumint + 1
            !IF(p_skip(dumint) .GE. 0) CYCLE
            If(p_d2(dumint,1) .LT. -1.0D7) CYCLE
            xp(k,1) = p_d2(dumint,1)
            xp(k,2) = p_d2(dumint,2)
            xp(k,3) = p_d2(dumint,3)

            vp(k,1) = p_d2(dumint,4)
            vp(k,2) = p_d2(dumint,5)
            vp(k,3) = p_d2(dumint,6)

            mp(k) = p_d2(dumint,7)
            ap(k) = p_d2(dumint,8)
            zp(k) = p_d2(dumint,9)

          ENDDO

          !!
          IF(matchok .GE. 0) dom_list2(g_num,i) = 1
          !!
          DEALLOCATE(p_i2, p_d2, p_skip)
        ENDDO

        DEALLOCATE(p_d, p_i)
      ENDDO
      !$OMP END PARALLEL DO

      js_order = 0

      !!!!!-- Matching Rate
      !$OMP PARALLEL DO default(shared) private(ind1, ind2, n0, n1, j) schedule(dynamic)
      DO i=1, n_gal
        ind1 = ind_b(i,1) + 1
        ind2 = ind_b(i,2) + 1

        n1 = 0
        n0 = ind2 - ind1 + 1
        DO j=ind1, ind2
          IF(mp(j) .GT. -1.0d7) n1 = n1 + 1
          IF(id(j) .LT. noptcl) n0 = n0 - 1
        ENDDO

        ind1 = ind_u(i,1) + 1
        ind2 = ind_u(i,2) + 1
        n0 = n0 + ind2 - ind1 + 1
        DO j=ind1, ind2
          IF(mp(j) .GT. -1.0d7) n1 = n1 + 1
          IF(id(j) .LT. noptcl) n0 = n0 - 1
        ENDDO

        rate(i) = float(n1) / float(n0)
      ENDDO
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO default(shared) schedule(static)
      DO i=1, n_gal
        DO j=1, n_mpi
          dom_list(i,j) = dom_list2(i,j)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      RETURN
      End

!!!!!
!! PARTICLE SORTING
!!!!!
      SUBROUTINE SORT_PTCL(p_dbl, p_int, nptcl)

      Implicit none
      Integer(kind=4) nptcl, i, j, ind, nptcl_dum0, nptcl_dum1

      Real(kind=8) p_dbl(nptcl,9), dumdbl(nptcl)
      Integer(kind=8) p_int(nptcl,3), dumid, p_int2(nptcl,2)
      Integer(kind=8) dumint(nptcl)

      INTEGER(KIND=4) first, last

      Do i=1, nptcl
        p_int2(i,1) = p_int(i,3)
        p_int2(i,2) = i
      Enddo

      nptcl_dum0 = 1
      nptcl_dum1 = nptcl
      Call quicksort(p_int2, nptcl, nptcl_dum0, nptcl_dum1)

      DO j=1,3
        DO i=1, nptcl
          dumint(i) = p_int(i,j)
        ENDDO

        DO i=1, nptcl
          p_int(i,j) = dumint(p_int2(i,2))
        ENDDO
      ENDDO

      DO j=1,9
        DO i=1, nptcl
          dumdbl(i) = p_dbl(i,j)
        ENDDO

        DO i=1, nptcl
          p_dbl(i,j) = dumdbl(p_int2(i,2))
        ENDDO
      ENDDO

      RETURN
      END

      !SUBROUTINE SORT_PTCL2(p_int, nptcl)

      !Implicit none
      !Integer(kind=4) nptcl, i, j, ind

      !Integer(kind=8) p_int(nptcl,3), dumid
      !Integer(kind=4) dumint

      !dumid = p_int(1,3)
      !DO i=1, nptcl
      !  dumid = p_int(i,3)
      !  ind = i

      !  DO j=i, nptcl
      !    IF(p_int(j,3) .le. dumid) THEN
      !      dumid = p_int(j,3)
      !      ind = j
      !    ENDIF
      !  ENDDO

      !  DO j=1, 3
      !    dumid = p_int(ind,j)
      !    p_int(ind,j) = p_int(i,j)
      !    p_int(i,j) = dumid
      !  ENDDO
      !ENDDO

      !RETURN
      !END


!!!!!
!! MATCHING
!!!!!
      SUBROUTINE match_hash(raw_int, raw_dbl, &
              gal_id, gal_dbl, &
              n_raw, n_gal, matchok, noptcl)
      IMPLICIT NONE
      INTEGER(KIND=4) n_raw, n_gal, matchok

      REAL(KIND=8) raw_dbl(n_raw,9)
      INTEGER(KIND=8) raw_int(n_raw, 3), raw_id(n_raw)
      INTEGER(KIND=8) noptcl

      REAL(KIND=8) gal_dbl(n_gal,9)
      INTEGER(KIND=8) gal_id(n_gal,1)
      !!-----
      !! LOCAL VARIABLEs
      !!-----
      INTEGER(KIND=4) i, j, k, l
      INTEGER(KIND=4) hash(n_raw), hash_next(n_raw)
      INTEGER(KIND=4) ind, i0

      DO i=1, n_raw
        raw_id(i) = raw_int(i,3)
      ENDDO

      hash      = -1
      hash_next = -1
      !!-----
      !! MAKE HASH TABLE
      !!-----
      DO i=1, n_raw
        ind = MOD(ABS(raw_id(i)), n_raw) + 1
        IF(hash(ind) .LT. 0) THEN
          hash(ind) = i
        ELSE
          i0 = hash(ind)
          DO WHILE (1 .EQ. 1)
            IF(hash_next(i0) .LE. 0) THEN
              hash_next(i0) = i
              EXIT
            ELSE
              i0 = hash_next(i0)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      !!-----
      !! MATCHING
      !!-----
      DO i=1, n_gal
        !PRINT *, i, n_raw, n_gal
        IF(gal_id(i,1) .LE. noptcl) CYCLE
        ind = MOD(ABS(gal_id(i,1)), n_raw) + 1
        i0  = hash(ind)
        IF(i0 .LE. 0) CYCLE
        DO WHILE(1 .EQ. 1)
          IF(gal_id(i,1) .EQ. raw_id(i0)) THEN
            DO j=1,9
              gal_dbl(i,j) = raw_dbl(i0,j)
            ENDDO
            matchok = 1
            EXIT
          ELSE
            IF(hash_next(i0) .GT. 0) THEN
              i0 = hash_next(i0)
            ELSE
              EXIT
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

      Subroutine match(raw_int, raw_dbl, &
              gal_id, gal_dbl, &
              n_raw, n_gal, matchok, larr, darr)

      implicit none
      integer(kind=4) n_raw, n_gal, matchok

      real(kind=8) raw_dbl(n_raw,9)
      integer(kind=8) raw_int(n_raw,3), raw_id(n_raw)

      real(kind=8) gal_dbl(n_gal,9)
      integer(kind=8) gal_id(n_gal, 1)

      integer(kind=4) larr(20)
      real(kind=8) darr(20)

!!!!!! Local Variables
      integer(kind=4) i, j, k, l, di
      integer(kind=4) ind0, ind1, mi, tid
      integer(kind=4) n_thread
      integer(kind=4) n_rep, n_test1, n_test2
      integer(kind=8), dimension(:), allocatable :: ind, id_rep
      !integer(kind=4) ind(larr(5))
      !integer(kind=8) id_rep(larr(5))
      logical gal_ch
        tid = matchok
      n_test1 = matchok
      n_test2 = 0
      matchok = -1
      !!!!
      DO i=1, n_raw
        raw_id(i) = raw_int(i,3)
      ENDDO

      DO i=1, n_gal
        IF(gal_id(i,1) .LE. raw_id(n_raw) .AND. &
                gal_id(i,1) .GE. raw_id(1)) THEN

          ind0 = 1
          ind1 = n_raw
          DO WHILE(ind1 - ind0 .GT. 10)
            l = int((ind1 + ind0)/2)

            IF(gal_id(i,1) .GT. raw_id(l)) ind0 = l
            IF(gal_id(i,1) .LT. raw_id(l)) ind1 = l
            IF(gal_id(i,1) .EQ. raw_id(l)) THEN
             ind0 = l; ind1 = l; mi = l
            ENDIF
          ENDDO

          IF(ind0 .NE. ind1) THEN
            mi = -1
            DO l=ind0, ind1
              IF(gal_id(i,1) .EQ. raw_id(l)) mi = l
            ENDDO
          ELSE
            mi = ind0
          ENDIF

          IF(mi .GE. 1) THEN
            DO j=1,9
              gal_dbl(i,j) = raw_dbl(mi,j)
            ENDDO
            matchok = 1
            n_test2 = n_test2 + 1
          ENDIF
        ENDIF

        !IF(gal_id(i,1) .EQ. 37827512) THEN
        !  PRINT *, mi, ind0, ind1
        !  STOP
        !ENDIF
      ENDDO

      !IF(n_test2 .GT. 0) PRINT *, n_test1, ' / ', n_test2
      RETURN

      !!!!

      !If(n_raw .le. 200) Then
      !  Do i=1, n_gal
      !    DO l=1, n_raw
      !      gal_ch = (raw_id(l) .eq. gal_id(i))
      !      If (gal_ch) then
      !        Do j=1,9
      !          gal_dbl(i,j) = raw_dbl(l,j)
      !        Enddo
      !      Endif
      !    ENDDO
      !  Enddo
      !Else
      !  n_rep = int(n_raw / 100)

      !  ALLOCATE(ind(1:n_rep))
      !  ALLOCATE(id_rep(1:n_rep))

      !  k = int(n_raw/(n_rep-1))
      !  Do i=2, n_rep-1
      !    ind(i) = k*i
      !  Enddo
      !  ind(1) = 1; ind(n_rep) = n_raw

      !  DO i=1, n_rep
      !    id_rep(i) = raw_id(ind(i))
      !  ENDDO
      !  id_rep(n_rep) = id_rep(n_rep) + 1


      !  DO i=1, n_gal
      !    di = -1
      !    Do j=1, n_rep-1
      !      IF(gal_id(i) .ge. id_rep(j) .and. &
      !          gal_id(i) .lt. id_rep(j+1)) Then
      !        di = j
      !        Exit
      !      Endif
      !    ENDDO

      !    IF(di .ge. 0) then
      !      DO l=ind(di), ind(di+1)
      !        gal_ch = (raw_id(l) .eq. gal_id(i))
      !        if (gal_ch) then
      !          Do j=1,9
      !            gal_dbl(i,j) = raw_dbl(l,j)
      !          Enddo
      !        endif
      !      ENDDO
      !    ENDIF
      !  ENDDO

      !  DEALLOCATE(ind, id_rep)
      !Endif

      !Return
      End subroutine

!!!!!
!! GET TARGETTED PARTICLE
!!!!!
      SUBROUTINE GET_PTCL(p_dbl, p_int, p_dbl2, p_int2, &
                      nbody, nptcl, ptype, dmp_mass)

      Implicit none
      Integer(kind=4) nbody, nptcl
      Real(kind=8) ptype, dmp_mass

      Real(kind=8) p_dbl(nptcl,9), p_dbl2(nbody,9)
      Integer(kind=8) p_int(nptcl,3), p_int2(nbody,3)

      Integer(kind=4) i, j, nn

      nn = 1
      Do i=1, nbody
        IF(ptype .gt. 0 .and. p_int2(i,2) .eq. 2) THEN
          DO j=1, 9
            p_dbl(nn,j) = p_dbl2(i,j)
          ENDDO

          DO j=1,3
            p_int(nn,j) = p_int2(i,j)
          ENDDO
          nn = nn + 1
        ENDIF
        IF(ptype .lt. 0 .and. abs(p_dbl2(i,7) - dmp_mass)/dmp_mass .lt. 1e-5) THEN
          DO j=1, 9
            p_dbl(nn,j) = p_dbl2(i,j)
          ENDDO

          DO j=1,3
            p_int(nn,j) = p_int2(i,j)
          ENDDO
          nn = nn + 1
        ENDIF
      ENDDO

      RETURN
      END


!!!!!
!! GET TARGETTED PARTICLE_YZiCS
!!!!!
      SUBROUTINE GET_PTCL_YZiCS(p_dbl, p_int, p_dbl2, p_int2, &
                      nbody, nptcl, ptype, dmp_mass)

      Implicit none
      Integer(kind=4) nbody, nptcl
      Real(kind=8) ptype, dmp_mass

      Real(kind=8) p_dbl(nptcl,9), p_dbl2(nbody,9)
      Integer(kind=8) p_int(nptcl,3), p_int2(nbody,3)

      Integer(kind=4) i, j, nn

      nn = 1
      Do i=1, nbody
        IF(ptype .gt. 0 .and. p_dbl2(i,8) .LT. 0) THEN
          DO j=1, 9
            p_dbl(nn,j) = p_dbl2(i,j)
          ENDDO

          DO j=1,3
            p_int(nn,j) = p_int2(i,j)
          ENDDO
          nn = nn + 1
        ENDIF
        IF(ptype .lt. 0 .and. abs(p_dbl2(i,7) - dmp_mass)/dmp_mass .lt. 1e-5) THEN
          DO j=1, 9
            p_dbl(nn,j) = p_dbl2(i,j)
          ENDDO

          DO j=1,3
            p_int(nn,j) = p_int2(i,j)
          ENDDO
          nn = nn + 1
        ENDIF
      ENDDO

      RETURN
      END

!!!!!
!! GET TARGETTED PARTICLE NUMBER
!!!!!
      SUBROUTINE GET_PTCL_NUM(p_m, p_fam, nbody, nptcl, ptype, dmp_mass)

      Implicit none
      Integer(kind=4) nptcl, nbody, i
      Real(kind=8) p_m(nbody), ptype, dmp_mass
      Integer(kind=8) p_fam(nbody)

      nptcl = 0
      Do i=1, nbody
        IF(ptype .gt. 0 .and. p_fam(i) .eq. 2) nptcl = nptcl + 1
        IF(ptype .lt. 0 .and. abs(p_m(i) - dmp_mass)/dmp_mass .lt. 1e-5) &
                nptcl = nptcl + 1
      ENDDO

      RETURN
      END

!!!!!
!! GET TARGETTED PARTICLE NUMBER_YZiCS
!!!!!
      SUBROUTINE GET_PTCL_NUM_YZiCS(p_m, p_age, nbody, nptcl, ptype, dmp_mass)

      Implicit none
      Integer(kind=4) nptcl, nbody, i
      Real(kind=8) p_m(nbody), p_age(nbody), ptype, dmp_mass

      nptcl = 0
      Do i=1, nbody
        IF(ptype .LT. 0 .AND. abs(p_m(i) - dmp_mass)/dmp_mass .lt. 1e-5) &
                nptcl = nptcl + 1

        IF(ptype .GT. 0 .AND. p_age(i) .LT. 0) &
                nptcl = nptcl + 1
      ENDDO

      RETURN
      END

!!!!!
!! READ PARTICLE
!!!!!
      SUBROUTINE RD_PART(dir_raw, n_snap, icpu, nbody, threadnum, &
               larr, p_dbl, p_int)

      Implicit none
      Integer(kind=4) larr(20)
      Integer(kind=4) icpu, n_snap, nbody, threadnum
      Character(LEN=larr(11)) dir_raw

      Real(kind=8) p_dbl(nbody,9), dumdbl
      Integer(kind=8) p_int(nbody,3)

      !! POS(3), VEL(3), MASS, AGE, METAL
      !! TAG, FAMILY, ID
!!!!! Local variables
      Integer(kind=4) i, j, k, uout
      Integer(kind=4) longint

      Real(kind=8) dum_dbl(nbody)
      Integer(kind=8) dum_int_ll(nbody)
      Integer(kind=4) dum_int(nbody)
      Integer(kind=1) dum_int_byte(nbody)

      Character*(100) fname, snap, domnum

      uout      = 100 + threadnum
      longint   = larr(20)

      write(snap, '(I5.5)') n_snap
      write(domnum, '(I5.5)') icpu
      fname = TRIM(dir_raw)//'output_'//TRIM(snap)//'/part_'//&
        TRIM(snap)//'.out'//TRIM(domnum)

      open(unit=uout, file=TRIM(fname), &
        form='unformatted', status='old')
      read(uout); read(uout); read(uout); read(uout)
      read(uout); read(uout); read(uout); read(uout)

      Do j=1, 3                         !! Position
        read(uout) dum_dbl
        Do i=1, nbody
          p_dbl(i,j) = dum_dbl(i)
        ENDDO
      ENDDO

      Do j=1, 3                         !! Velocity
        read(uout) dum_dbl
        Do i=1, nbody
          p_dbl(i,j+3) = dum_dbl(i)
        ENDDO
      ENDDO

      read(uout) dum_dbl                !! Mass
      Do i=1, nbody
        p_dbl(i,7) = dum_dbl(i)
      ENDDO


      IF(longint .le. 10) THEN
        read(uout) dum_int_ll
        Do i=1, nbody
          p_int(i,3) = dum_int_ll(i)
        ENDDO
      ELSE
        read(uout) dum_int
        Do i=1, nbody
          p_int(i,3) = dum_int(i)
        ENDDO
      ENDIF

      read(uout) dum_int

      read(uout) dum_int_byte
      Do i=1, nbody
        p_int(i,2) = dum_int_byte(i)
        IF(p_int(i,2) .gt. 100) p_int(i,2) = p_int(i,2) - 255
      ENDDO

      read(uout) dum_int_byte
      Do i=1, nbody
        p_int(i,1) = dum_int_byte(i)
        IF(p_int(i,1) .gt. 100) p_int(i,1) = p_int(i,1) - 255
      ENDDO

      read(uout) dum_dbl                !! Age
      Do i=1, nbody
        p_dbl(i,8) = dum_dbl(i)
      ENDDO

      read(uout) dum_dbl                !! Metallicity
      Do i=1, nbody
        p_dbl(i,9) = dum_dbl(i)
      ENDDO

      close(uout)
      RETURN
      END

!!!!!
!! READ PARTICLE_YZICS
!!!!!
      SUBROUTINE RD_PART_YZiCS(dir_raw, n_snap, icpu, nbody, threadnum,&
        larr, p_dbl, p_int)

      Implicit none
      Integer(kind=4) larr(20)
      Integer(kind=4) icpu, nbody, threadnum, n_snap
      Character(LEN=larr(11)) dir_raw

      Real(kind=8) p_dbl(nbody,9), dumdbl
      Integer(kind=8) p_int(nbody,3)

      !! POS(3), VEL(3), MASS, AGE, METAL
      !! TAG, FAMILY, ID
!!!!! Local variables
      Integer(kind=4) i, j, k, uout
      Integer(kind=4) longint

      Real(kind=8) dum_dbl(nbody)
      Integer(kind=8) dum_int_ll(nbody)
      Integer(kind=4) dum_int(nbody)
      Integer(kind=1) dum_int_byte(nbody)

      Character*(100) fname, snap, domnum

      longint   = larr(20)
      uout      = 100 + threadnum

      write(snap, '(I5.5)') n_snap
      write(domnum, '(I5.5)') icpu
      fname = TRIM(dir_raw)//'output_'//TRIM(snap)//'/part_'//&
        TRIM(snap)//'.out'//TRIM(domnum)

      open(unit=uout, file=TRIM(fname), &
        form='unformatted', status='old')
      read(uout); read(uout); read(uout); read(uout)
      read(uout); read(uout); read(uout); read(uout)

      Do j=1, 3                         !! Position
        read(uout) dum_dbl
        Do i=1, nbody
          p_dbl(i,j) = dum_dbl(i)
        ENDDO
      ENDDO

      Do j=1, 3                         !! Velocity
        read(uout) dum_dbl
        Do i=1, nbody
          p_dbl(i,j+3) = dum_dbl(i)
        ENDDO
      ENDDO

      read(uout) dum_dbl                !! Mass
      Do i=1, nbody
        p_dbl(i,7) = dum_dbl(i)
      ENDDO


      IF(longint .le. 10) THEN          !! ID
        read(uout) dum_int_ll
        Do i=1, nbody
          p_int(i,3) = dum_int_ll(i)
        ENDDO
      ELSE
        read(uout) dum_int
        Do i=1, nbody
          p_int(i,3) = dum_int(i)
        ENDDO
      ENDIF

      read(uout) dum_int                !! LEVEL

      !read(uout) dum_int_byte
      !Do i=1, nbody
      !  p_int(i,2) = dum_int_byte(i)
      !  IF(p_int(i,2) .gt. 100) p_int(i,2) = p_int(i,2) - 255
      !ENDDO

      !read(uout) dum_int_byte
      !Do i=1, nbody
      !  p_int(i,1) = dum_int_byte(i)
      !  IF(p_int(i,1) .gt. 100) p_int(i,1) = p_int(i,1) - 255
      !ENDDO

      read(uout) dum_dbl                !! Age
      Do i=1, nbody
        p_dbl(i,8) = dum_dbl(i)
      ENDDO

      read(uout) dum_dbl                !! Metallicity
      Do i=1, nbody
        p_dbl(i,9) = dum_dbl(i)
      ENDDO

      close(uout)
      RETURN
      END
!!!!!
!! GET NBODY
!!!!!
      SUBROUTINE RD_PART_NBODY(dir_raw, n_snap, icpu, &
                threadnum, nbody, larr)

      Implicit none
      Integer(kind=4) larr(20)
      Integer(kind=4) n_snap, threadnum, icpu
      Character(LEN=larr(11)) dir_raw

!!!!! Local variables
      Character*(100) fname, snap, domnum
      Integer(kind=4) uout, nbody
      uout      = 100 + threadnum

      write(snap, '(I5.5)') n_snap
      write(domnum, '(I5.5)') icpu
      fname = TRIM(dir_raw)//'output_'//TRIM(snap)//'/part_'//&
        TRIM(snap)//'.out'//TRIM(domnum)

      open(unit=uout, file=TRIM(fname), form='unformatted', status='old')
      read(uout)
      read(uout)
      read(uout) nbody
      close(uout)

      RETURN
      END
!!!!!
!! GAL LIST
!!!!!
      SUBROUTINE GAL_LIST(dom_list, n_gal, n_mpi, mpi_num, gals, g_num)

      Implicit none
      Integer(kind=4) n_gal, n_mpi
      Integer(kind=4) dom_list(n_gal, n_mpi)
      Integer(kind=4) mpi_num
      Integer(kind=4) gals(n_gal)
      Integer(kind=4) g_num
      Integer(kind=4) i, j, k

      Do i=1, n_gal
        gals(i) = -1
      ENDDO

      g_num = -1
      Do i=1, n_gal
        IF(DOM_LIST(i,mpi_num) .ge. 0) THEN
          IF(g_num .lt. 0) g_num = 1
          gals(g_num) = i
          g_num = g_num + 1
        ENDIF
      ENDDO

      IF(g_num .ge. 0) g_num = g_num - 1
      RETURN
      END

!!!!!
!! QUICK SORT
!!!!!
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
      recursive subroutine quicksort(a, nn, first, last)
        implicit none
        integer*8  a(nn,2), x, t, n
        integer first, last
        integer i, j, nn
      
        x = a( (first+last) / 2 , 1)
        i = first
        j = last
        do
           do while (a(i,1) < x)
              i=i+1
           end do
           do while (x < a(j,1))
              j=j-1
           end do
           IF(i .GE. j) exit
           t = a(i,1);  a(i,1) = a(j,1);  a(j,1) = t
           n = a(i,2);  a(i,2) = a(j,2);  a(j,2) = n
           i=i+1
           j=j-1
        end do
        if (first < i-1) call quicksort(a, nn, first, i-1)
        if (j+1 < last)  call quicksort(a, nn, j+1, last)
        return
      end subroutine quicksort
