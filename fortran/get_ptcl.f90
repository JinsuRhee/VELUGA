!234567
      Subroutine get_ptcl(larr, darr, dir_raw, &
                ID, ptcl, dom_list)

      USE omp_lib

      Implicit none
      Integer(kind=4) larr(20)
      Real(kind=8) darr(20)

      Integer(kind=8) ID(larr(1))
      Real(kind=8) ptcl(larr(1),9)
      Integer(kind=4) dom_list(larr(2))
      Character*(larr(11)) dir_raw

!!!!!! Local Variables
      Integer(kind=4) n_thread, n_ptcl, n_snap, n_raw, n_dom, n_str
      INTEGER(KIND=4) longint, n_star
      Integer(kind=4) i, j, k, l, ind0, ind1, mi, horg
      REAL(KIND=8), allocatable, dimension(:,:) :: raw_dbl, raw_dbl2
      INTEGER(KIND=8), allocatable, dimension(:,:) :: raw_int, raw_int2
      INTEGER(KIND=4) dom_index(larr(2)), p_index(larr(4))
      REAL(KIND=8) dmp_mass

      n_ptcl    = larr(1)
      n_dom     = larr(2)
      n_snap    = larr(3)
      n_thread  = larr(4)
      n_str     = larr(11)
      longint   = larr(20)
      horg      = larr(12)
      IF(horg .EQ.0) horg = 1   !! star particle search as a default
      dmp_mass  = darr(12)

      CALL OMP_SET_NUM_THREADS(n_thread)

!!!!! RD PART

      CALL RD_PART_NBODY(dir_raw, dom_list, n_dom, n_raw, n_str, n_snap, &
             dom_index, n_thread)

      ALLOCATE(raw_dbl(1:n_raw,1:9))
      ALLOCATE(raw_int(1:n_raw,1:2))

      IF(larr(19) .GT. 10) &
        CALL RD_PART(dir_raw, dom_list, n_dom, n_raw, n_str, n_snap, &
              raw_dbl, raw_int, longint, dom_index, n_thread)

      IF(larr(19) .LT. 10) &
        CALL RD_PART_YZiCS(dir_raw, dom_list, n_dom, n_raw, n_str, n_snap, &
              raw_dbl, raw_int, longint, dom_index, n_thread)

      IF(larr(19) .GT. 10) &
        CALL GET_PTCL_NUM(raw_int, raw_dbl, n_raw, n_star, n_thread, &
                horg, dmp_mass, p_index)
      IF(larr(19) .LT. 10) &
        CALL GET_PTCL_NUM_YZiCS(raw_dbl, n_raw, n_star, n_thread, &
                horg, dmp_mass, p_index)

      ALLOCATE(raw_dbl2(1:n_star,1:9))
      ALLOCATE(raw_int2(1:n_star,1:1))

      IF(larr(19) .GT. 10) &
        CALL GET_STAR_PTCL(raw_dbl, raw_int, raw_dbl2, raw_int2, &
              n_raw, n_star, horg, dmp_mass, n_thread, p_index)
      IF(larr(19) .LT. 10) &
        CALL GET_STAR_PTCL_YZiCS(raw_dbl, raw_int, raw_dbl2, raw_int2, &
              n_raw, n_star, horg, dmp_mass, n_thread, p_index)

      DEALLOCATE(raw_dbl, raw_int)

      IF(larr(18) .GE. 10) THEN !Return for calling raw particles
        !$OMP PARALLEL DO default(shared)
        DO i=1, n_star
          DO j=1, 9
            ptcl(i,j) = raw_dbl2(i,j)
          ENDDO
          id(i) = raw_int2(i,1)
        ENDDO
        !$OMP END PARALLEL DO
        larr(17) = n_star
        RETURN
      ENDIF

      CALL SORT_PTCL(raw_dbl2, raw_int2, n_star)

      !!-----
      !! MATCHING
      !!  *) To do list: test w/ the hash search
      !!-----
      !$OMP PARALLEL DO default(shared) private(ind0, ind1, l, mi) schedule(dynamic)
      DO i=1, n_ptcl
        IF(ID(i) .LE. raw_int2(n_star,1) .AND. &
                ID(i) .GE. raw_int2(1,1)) THEN

          ind0 = 1
          ind1 = n_star
          DO WHILE (ind1 - ind0 .GT. 10)
            l = int((ind1 + ind0) / 2)
            IF(ID(i) .GT. raw_int2(l,1)) ind0 = l
            IF(ID(i) .LT. raw_int2(l,1)) ind1 = l
            IF(ID(i) .EQ. raw_int2(l,1)) THEN
              ind0 = l; ind1 = l; mi = l
            ENDIF
          ENDDO

          IF(ind0 .NE. ind1) THEN
            mi = -1
            DO l=ind0, ind1
              IF(ID(i) .EQ. raw_int2(l,1)) mi = l
            ENDDO
          ELSE
            mi = ind0
          ENDIF


          IF(mi .GE. 1) THEN
            DO j=1,9
              ptcl(i,j) = raw_dbl2(mi,j)
            ENDDO
          ENDIF 
        ENDIF
      ENDDO
      !$OMP END PARALLEL DO


      RETURN
      END

!!!!!
!! SORT BY PARTICLE ID
!!!!!
      SUBROUTINE SORT_PTCL(p_dbl, p_int, nptcl)

      Implicit none
      Integer(kind=4) nptcl, i, j, ind, nptcl_dum0, nptcl_dum1

      Real(kind=8) p_dbl(nptcl,9), dumdbl(nptcl)
      Integer(kind=8) p_int(nptcl,1), dumid, p_int2(nptcl,2)
      Integer(kind=8) dumint(nptcl)

      INTEGER(KIND=4) first, last

      Do i=1, nptcl
        p_int2(i,1) = p_int(i,1)
        p_int2(i,2) = i
      Enddo

      nptcl_dum0 = 1
      nptcl_dum1 = nptcl
      Call quicksort(p_int2, nptcl, nptcl_dum0, nptcl_dum1)

      DO i=1, nptcl
        dumint(i) = p_int(i,1)
      ENDDO

      DO i=1, nptcl
        p_int(i,1) = dumint(p_int2(i,2))
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
!!!!!
!! GET TARGETTED PARTICLE
!!!!!
      SUBROUTINE GET_STAR_PTCL(p_dbl, p_int, p_dbl2, p_int2, &
                      n_raw, n_star, horg, dmp_mass, n_thread, p_index)

      Implicit none
      Integer(kind=4) n_raw, n_star, horg, n_thread, p_index(n_thread)

      Real(kind=8) p_dbl(n_raw,9), p_dbl2(n_star,9), dmp_mass
      Integer(kind=8) p_int(n_raw,2), p_int2(n_star,1)

      Integer(kind=4) i, j, nn, res, t_ind0(n_thread), t_ind

      res       = (n_raw*1.0) / (n_thread*1.0) + 1
      t_ind0    = 1
      IF(horg .GT. 0.) THEN !STAR
        !$OMP PARALLEL DO default(shared) schedule(static, res) &
        !$OMP & private(t_ind)
        Do i=1, n_raw
          IF(p_int(i,2) .eq. 2) THEN
            t_ind = INT((i-1)/res) + 1
            IF(t_ind .GT. n_thread) t_ind = n_thread

            DO j=1, 9
              p_dbl2(p_index(t_ind) + t_ind0(t_ind),j) = p_dbl(i,j)
            ENDDO
            p_int2(p_index(t_ind) + t_ind0(t_ind),1) = p_int(i,1)

            t_ind0(t_ind) = t_ind0(t_ind) + 1
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ELSE IF (horg .LT. 0) THEN !DM
        !$OMP PARALLEL DO default(shared) schedule(static, res) &
        !$OMP & private(t_ind)
        Do i=1, n_raw
          IF(p_int(i,2) .eq. 1) THEN
            IF(ABS(p_dbl(i,7)-dmp_mass)/dmp_mass .LT. 1e-5) THEN
              t_ind = INT((i-1)/res) + 1
              IF(t_ind .GT. n_thread) t_ind = n_thread

              DO j=1, 9
                p_dbl2(p_index(t_ind) + t_ind0(t_ind),j) = p_dbl(i,j)
              ENDDO
              p_int2(p_index(t_ind) + t_ind0(t_ind),1) = p_int(i,1)

              t_ind0(t_ind) = t_ind0(t_ind) + 1
            ENDIF
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF

      RETURN
      END
!!!!!
!! GET TARGETTED PARTICLE_YZiCS
!!!!!
      SUBROUTINE GET_STAR_PTCL_YZiCS(p_dbl, p_int, p_dbl2, p_int2, &
                      n_raw, n_star, horg, dmp_mass, n_thread, p_index)

      USE omp_lib
      Implicit none
      Integer(kind=4) n_raw, n_star, horg, n_thread, p_index(n_thread)

      Real(kind=8) p_dbl(n_raw,9), p_dbl2(n_star,9), dmp_mass
      Integer(kind=8) p_int(n_raw,2), p_int2(n_star,1)

      Integer(kind=4) i, j, nn, res, t_ind, t_num, t_ind0(n_thread)

      res       = (n_raw*1.0) / (n_thread*1.0) + 1
      t_ind0 = 1
      IF(horg .GT. 0) THEN !STAR
        !$OMP PARALLEL DO default(shared) schedule(static, res) &
        !$OMP & private(t_ind)
        Do i=1, n_raw
          IF(p_dbl(i,8) .LT. 0) THEN
            t_ind = INT((i-1)/res) + 1
            IF(t_ind .GT. n_thread) t_ind = n_thread

            DO j=1, 9
              p_dbl2(p_index(t_ind) + t_ind0(t_ind),j) = p_dbl(i,j)
            ENDDO
            p_int2(p_index(t_ind) + t_ind0(t_ind),1) = p_int(i,1)

            t_ind0(t_ind) = t_ind0(t_ind) + 1
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ELSE IF (horg .LT. 0) THEN !DM
        !$OMP PARALLEL DO default(shared) schedule(static, res) &
        !$OMP & private(t_ind)
        Do i=1, n_raw
          IF(p_dbl(i,8) .EQ. 0) THEN
            t_ind = INT((i-1)/res) + 1
            IF(t_ind .GT. n_thread) t_ind = n_thread

            IF(ABS(p_dbl(i,7)-dmp_mass)/dmp_mass .LT. 1e-5) THEN
              DO j=1, 9
                p_dbl2(p_index(t_ind) + t_ind0(t_ind),j) = p_dbl(i,j)
              ENDDO
              p_int2(p_index(t_ind) + t_ind0(t_ind),1) = p_int(i,1)

              t_ind0(t_ind) = t_ind0(t_ind) + 1
            ENDIF
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
      RETURN
      END
!!!!!
!! GET PARTICLE NUM
!!!!!
      SUBROUTINE GET_PTCL_NUM(rint, rdbl, n_raw, n_star, n_thread, &
                  horg, dmp_mass, p_ind)

      USE omp_lib
      IMPLICIT NONE
      INTEGER(KIND=4) n_raw, n_star, n_thread, p_ind(n_thread)
      INTEGER(KIND=8) rint(n_raw, 2)
      REAL(KIND=8) rdbl(n_raw,9)

      INTEGER(KIND=4) i, j, k, horg, t_ind, parts(n_thread)
      INTEGER(KIND=4) res
      REAL(KIND=8) dmp_mass

      CALL OMP_SET_NUM_THREADS(n_thread)
      n_star = 0
      parts  = 0
      res    = (n_raw*1.0) / (n_thread*1.0) + 1
      IF(horg .GT. 0.) THEN ! STAR
        !$OMP PARALLEL DO default(shared) schedule(static, res) & 
        !$OMP private(t_ind) reduction(+:n_star, parts)
        DO i=1, n_raw
          t_ind = INT((i-1)/res) + 1
          IF(t_ind .GT. n_thread) t_ind = n_thread

          IF(rint(i,2) .EQ. 2) THEN
            n_star = n_star + 1
            parts(t_ind) = parts(t_ind) + 1
          ENDIF

        ENDDO
        !$OMP END PARALLEL DO

      ELSE IF (horg .LT. 0.) THEN !DM
        !$OMP PARALLEL DO default(shared) schedule(static, res) &
        !$OMP & private(t_ind) reduction(+:n_star, parts)
        DO i=1, n_raw
          t_ind = INT((i-1)/res) + 1
          IF(t_ind .GT. n_thread) t_ind = n_thread

          IF(rint(i,2) .EQ. 1) THEN
            IF(ABS(rdbl(i,7) - dmp_mass)/dmp_mass .LT. 1e-5) THEN
                n_star = n_star + 1
                parts(t_ind) = parts(t_ind) + 1
            ENDIF
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF

      p_ind     = 0
      IF(n_thread .GT. 1) THEN
        DO i=2, n_thread
          DO j=1, i-1
            p_ind(i) = p_ind(i) + parts(j)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END

!!!!!
!! GET PARTICLE NUM_YZiCS
!!!!!
      SUBROUTINE GET_PTCL_NUM_YZiCS(age, n_raw, n_star, n_thread, &
                      horg, dmp_mass, p_ind)

      USE omp_lib
      IMPLICIT NONE
      INTEGER(KIND=4) n_raw, n_star, n_thread, p_ind(n_thread)
      REAL(KIND=8) age(n_raw,9), dmp_mass

      INTEGER(KIND=4) i, j, k, horg, t_ind, parts(n_thread), tt
      INTEGER(KIND=4) res

      CALL OMP_SET_NUM_THREADS(n_thread)
      n_star = 0
      parts  = 0
      res       = (n_raw*1.0) / (n_thread*1.0) + 1

      IF(horg .GT. 0) THEN !STAR
        !$OMP PARALLEL DO default(shared) schedule(static, res) &
        !$OMP & private(t_ind, tt) reduction(+:n_star, parts)
        DO i=1, n_raw
          t_ind = INT((i-1)/res) + 1
          IF(t_ind .GT. n_thread) t_ind = n_thread

          IF(age(i,8) .LT. 0) THEN
            n_star = n_star + 1
            parts(t_ind) = parts(t_ind) + 1
          ENDIF

        ENDDO
        !$OMP END PARALLEL DO
      ELSE IF (horg .LT. 0) THEN
        !$OMP PARALLEL DO default(shared) schedule(static, res) & 
        !$OMP & private(t_ind) reduction(+:n_star, parts)
        DO i=1, n_raw
          t_ind = INT((i-1)/res) + 1
          IF(t_ind .GT. n_thread) t_ind = n_thread

          IF(age(i,8) .EQ. 0) THEN
           IF(ABS(age(i,7) - dmp_mass)/dmp_mass .LT. 1e-5) THEN
             n_star = n_star + 1
             parts(t_ind) = parts(t_ind) + 1
           ENDIF
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF

      p_ind = 0
      IF(n_thread .GT. 1) THEN
        DO i=2, n_thread
          DO j=1, i-1
            p_ind(i) = p_ind(i) + parts(j)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END

!!!!!
!! RD_PART
!!!!!
      SUBROUTINE RD_PART(dir_raw, dom_list, &
                n_dom, n_raw, n_str, n_snap, &
                rdbl, rint, longint, dom_ind, n_thread)
      USE omp_lib
      IMPLICIT NONE
      INTEGER(KIND=4) n_dom, n_raw, n_str, n_snap, n_thread
      INTEGER(KIND=4) dom_list(n_dom), longint
      CHARACTER*(n_str) dir_raw
      REAL(KIND=8) rdbl(n_raw,9)
      INTEGER(KIND=8) rint(n_raw,2)
      INTEGER(KIND=4) dom_ind(n_dom)

!!!!! Local variables

      CHARACTER*(100) fname, snap, domnum
      INTEGER(KIND=4) uout, icpu, nbody, n_thread2
      INTEGER(KIND=4) i, j, k, nn
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: dum_dbl
      INTEGER(KIND=8), ALLOCATABLE, DIMENSION(:) :: dum_int_ll
      INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: dum_int
      INTEGER(KIND=1), ALLOCATABLE, DIMENSION(:) :: dum_int_byte

      n_thread2 = n_thread
      IF(n_thread .GT. n_dom) n_thread2 = n_dom
      CALL OMP_SET_NUM_THREADS(n_thread2)

      nn = 0
      !$OMP PARALLEL DO default(shared) &
      !$OMP & schedule(static, 10) &
      !$OMP & private(snap, domnum, fname, uout, dum_dbl, dum_int_ll) &
      !$OMP & private(dum_int, dum_int_byte, j, i, nbody)
      DO k=1, n_dom
        WRITE(snap, '(I5.5)') n_snap
        WRITE(domnum, '(I5.5)') dom_list(k)
        fname = TRIM(dir_raw)//'output_'//TRIM(snap)//'/part_'//&
          TRIM(snap)//'.out'//TRIM(domnum)

        uout    = OMP_GET_THREAD_NUM() 
        uout    = uout + 10

        open(newunit=uout, file=fname, form='unformatted', status='old')
        read(uout); read(uout); read(uout) nbody; read(uout)
        read(uout); read(uout); read(uout); read(uout)

        ALLOCATE(dum_dbl(1:nbody))
        ALLOCATE(dum_int_ll(1:nbody))
        ALLOCATE(dum_int(1:nbody))
        ALLOCATE(dum_int_byte(1:nbody))

        Do j=1, 3                         !! Position
          read(uout) dum_dbl
          Do i=1, nbody
            !rdbl(i+nn,j) = dum_dbl(i)
            rdbl(dom_ind(k)+i,j) = dum_dbl(i)
          ENDDO
        ENDDO

        Do j=1, 3                         !! Velocity
          read(uout) dum_dbl
          Do i=1, nbody
            !rdbl(i+nn,j+3) = dum_dbl(i)
            rdbl(dom_ind(k)+i,j+3) = dum_dbl(i)
          ENDDO
        ENDDO

        read(uout) dum_dbl                !! Mass
        Do i=1, nbody
          !rdbl(i+nn,7) = dum_dbl(i)
          rdbl(dom_ind(k)+i,7) = dum_dbl(i)
        ENDDO

        IF(longint .ge. 10) THEN
          read(uout) dum_int_ll
          Do i=1, nbody
            !rint(i+nn,1) = dum_int_ll(i)
            rint(dom_ind(k)+i,1) = dum_int_ll(i)
          ENDDO
        ELSE
          read(uout) dum_int
          Do i=1, nbody
            !rint(i+nn,1) = dum_int(i)
            rint(dom_ind(k)+i,1) = dum_int(i)
          ENDDO
        ENDIF

        read(uout) dum_int

        read(uout) dum_int_byte
        Do i=1, nbody
          rint(dom_ind(k)+i,2) = dum_int_byte(i)
          IF(rint(dom_ind(k)+i,2) .gt. 100) rint(dom_ind(k)+i,2) = rint(dom_ind(k)+i,2) - 255
          !rint(i+nn,2) = dum_int_byte(i)
          !IF(rint(i+nn,2) .gt. 100) rint(i+nn,2) = rint(i+nn,2) - 255
        ENDDO

        read(uout) dum_int_byte

        read(uout) dum_dbl                !! Age
        Do i=1, nbody
          !rdbl(i+nn,8) = dum_dbl(i)
          rdbl(dom_ind(k)+i,8) = dum_dbl(i)
        ENDDO

        read(uout) dum_dbl                !! Metallicity
        Do i=1, nbody
          !rdbl(i+nn,9) = dum_dbl(i)
          rdbl(dom_ind(k)+i,9) = dum_dbl(i)
        ENDDO

        DEALLOCATE(dum_dbl, dum_int_ll)
        DEALLOCATE(dum_int, dum_int_byte)
        !nn = nn + nbody
        CLOSE(uout)
      ENDDO
      !$OMP END PARALLEL DO

      RETURN
      END 

!!!!!
!! RD_PART_YZiCS
!!!!!
      SUBROUTINE RD_PART_YZiCS(dir_raw, dom_list, &
                n_dom, n_raw, n_str, n_snap, &
                rdbl, rint, longint, dom_ind, n_thread)
      USE omp_lib
      IMPLICIT NONE
      INTEGER(KIND=4) n_dom, n_raw, n_str, n_snap, n_thread
      INTEGER(KIND=4) dom_list(n_dom), longint
      CHARACTER*(n_str) dir_raw
      REAL(KIND=8) rdbl(n_raw,9)
      INTEGER(KIND=8) rint(n_raw,2)
      INTEGER(KIND=4) dom_ind(n_dom)

!!!!! Local variables

      CHARACTER*(100) fname, snap, domnum
      INTEGER(KIND=4) uout, icpu, nbody, n_thread2
      INTEGER(KIND=4) i, j, k, nn
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: dum_dbl
      INTEGER(KIND=8), ALLOCATABLE, DIMENSION(:) :: dum_int_ll
      INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: dum_int
      INTEGER(KIND=1), ALLOCATABLE, DIMENSION(:) :: dum_int_byte

      n_thread2 = n_thread
      IF(n_thread .GT. n_dom) n_thread2 = n_dom
      CALL OMP_SET_NUM_THREADS(n_thread2)

      nn = 0
      !$OMP PARALLEL DO default(shared) &
      !$OMP & schedule(static, 10) &
      !$OMP & private(snap, domnum, fname, uout, dum_dbl, dum_int_ll) &
      !$OMP & private(dum_int, dum_int_byte, j, i, nbody)
      DO k=1, n_dom
        WRITE(snap, '(I5.5)') n_snap
        WRITE(domnum, '(I5.5)') dom_list(k)
        fname = TRIM(dir_raw)//'output_'//TRIM(snap)//'/part_'//&
          TRIM(snap)//'.out'//TRIM(domnum)

        uout    = OMP_GET_THREAD_NUM()
        uout    = uout + 10

        open(newunit=uout, file=fname, form='unformatted', status='old')
        read(uout); read(uout); read(uout) nbody; read(uout)
        read(uout); read(uout); read(uout); read(uout)

        ALLOCATE(dum_dbl(1:nbody))
        ALLOCATE(dum_int_ll(1:nbody))
        ALLOCATE(dum_int(1:nbody))
        ALLOCATE(dum_int_byte(1:nbody))

        Do j=1, 3                         !! Position
          read(uout) dum_dbl
          Do i=1, nbody
            !rdbl(i+nn,j) = dum_dbl(i)
            rdbl(dom_ind(k)+i,j) = dum_dbl(i)
          ENDDO
        ENDDO

        Do j=1, 3                         !! Velocity
          read(uout) dum_dbl
          Do i=1, nbody
            !rdbl(i+nn,j+3) = dum_dbl(i)
            rdbl(dom_ind(k)+i,j+3) = dum_dbl(i)
          ENDDO
        ENDDO

        read(uout) dum_dbl                !! Mass
        Do i=1, nbody
          !rdbl(i+nn,7) = dum_dbl(i)
          rdbl(dom_ind(k)+i,7) = dum_dbl(i)
        ENDDO

        IF(longint .ge. 10) THEN          !! ID
          read(uout) dum_int_ll
          Do i=1, nbody
            !rint(i+nn,1) = dum_int_ll(i)
            rint(dom_ind(k)+i,1) = dum_int_ll(i)
          ENDDO
        ELSE
          read(uout) dum_int
          Do i=1, nbody
            !rint(i+nn,1) = dum_int(i)
            rint(dom_ind(k)+i,1) = dum_int(i)
          ENDDO
        ENDIF

        read(uout) dum_int                !! LEVEL

        !read(uout) dum_int_byte
        !Do i=1, nbody
        !  rint(i+nn,2) = dum_int_byte(i)
        !  IF(rint(i+nn,2) .gt. 100) rint(i+nn,2) = rint(i+nn,2) - 255
        !ENDDO

        !read(uout) dum_int_byte

        read(uout) dum_dbl                !! Age
        Do i=1, nbody
          !rdbl(i+nn,8) = dum_dbl(i)
          rdbl(dom_ind(k)+i,8) = dum_dbl(i)
        ENDDO

        read(uout) dum_dbl                !! Metallicity
        Do i=1, nbody
          !rdbl(i+nn,9) = dum_dbl(i)
          rdbl(dom_ind(k)+i,9) = dum_dbl(i)
        ENDDO

        DEALLOCATE(dum_dbl, dum_int_ll)
        DEALLOCATE(dum_int, dum_int_byte)
        !nn = nn + nbody
        CLOSE(uout)
      ENDDO
      !$OMP END PARALLEL DO

      RETURN
      END 

!!!!!
!! GET NBODY
!!!!!
      SUBROUTINE RD_PART_NBODY(dir_raw, dom_list, &
              n_dom, n_raw, n_str, n_snap, dom_ind, n_thread)
      USE omp_lib
      Implicit none
      Integer(kind=4) n_raw, n_dom, n_str, n_snap, dom_list(n_dom)
      INTEGER(KIND=4) n_thread, dom_ind(n_dom)
      Character*(n_str) dir_raw

!!!!! Local variables
      Character*(100) fname, snap, domnum
      Integer(kind=4) uout, icpu, nbody, parts(n_dom)
      INTEGER(kind=4) i, j

      CALL OMP_SET_NUM_THREADS(n_thread)
      n_raw = 0
      !$OMP PARALLEL DO default(shared) schedule(static) &
      !$OMP & private(snap, domnum, fname, nbody, uout) &
      !$OMP reduction(+:n_raw)
      DO i=1, n_dom
        write(snap, '(I5.5)') n_snap
        write(domnum, '(I5.5)') dom_list(i)
        fname = TRIM(dir_raw)//'output_'//TRIM(snap)//'/part_'//&
          TRIM(snap)//'.out'//TRIM(domnum)

        uout    = OMP_GET_THREAD_NUM()
        uout    = uout + 10
        open(newunit=uout, file=fname, form='unformatted', status='old')
        read(uout)
        read(uout)
        read(uout) nbody
        close(uout)

        n_raw = n_raw + nbody
        parts(i)      = nbody
      ENDDO
      !$OMP END PARALLEL DO

      dom_ind = 0
      IF(n_dom .GT. 1) THEN
        DO i=2, n_dom
          DO j=1, i-1
            dom_ind(i)      = dom_ind(i) + parts(j)
          ENDDO
        ENDDO
      ENDIF

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
