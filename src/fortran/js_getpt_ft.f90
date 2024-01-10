!!---------------------------------------------
!! MAIN
!!---------------------------------------------

!234567
      SUBROUTINE js_getpt_ft(larr, darr, pos, mm, pot, force)

      USE omp_lib
      USE js_kdtree

      IMPLICIT NONE

      REAL(KIND=8) darr(20)
      INTEGER(KIND=4) larr(20)

      REAL(KIND=8) pos(larr(1),larr(2))
      REAL(KIND=8) mm(larr(1)), pot(larr(1)), force(larr(1))

      TYPE(nodetype), DIMENSION(:), ALLOCATABLE :: root
      !!-----
      !! LOCAL VARIABLES
      !!-----
      INTEGER(KIND=4) i, j, k, l, m

      INTEGER(KIND=4) n_ptcl, n_thread, n_dim, p_type, n_leaf
      INTEGER(KIND=4) d_type, v_type, bsize, e_type
      INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: orgind
      INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: recind
      REAL(KIND=8) Gconst, dummy, dx, time(10)
      REAL(KIND=8) dummy_v(8), bnd(8,3)
      REAL(KIND=8) dummy_vf(8)
      INTEGER(KIND=4) bs, be
      TYPE(nodetype), DIMENSION(:), ALLOCATABLE :: lf
      !TYPE(dat), DIMENSION(:), ALLOCATABLE :: part
      TYPE(infotype) :: tree_set

      n_ptcl    = larr(1)
      n_dim     = larr(2)
      n_thread  = larr(3)
      p_type    = larr(4)
      e_type    = larr(5)

      d_type    = larr(11)      ! DIMENSION TYPE
      v_type    = larr(12)      ! VALUE TYP
      bsize     = larr(13)
      Gconst    = darr(1)

      CALL OMP_SET_NUM_THREADS(n_thread)

      !!-----
      !! GET TREE
      !!-----
      IF(ALLOCATED(orgind)) DEALLOCATE(orgind)
      IF(ALLOCATED(recind)) DEALLOCATE(recind)
      ALLOCATE(orgind(1:n_ptcl))
      ALLOCATE(recind(1:n_ptcl))
      DO i=1, n_ptcl
        orgind(i) = i
        recind(i) = i
      ENDDO

      !IF(ALLOCATED(part)) DEALLOCATE(part)
      !ALLOCATE(part(1:n_ptcl))
      !DO i=1, n_dim
      !  part%pos(i) = pos(:,i)
      !ENDDO
      !part%mm = mm(:)

      time(1)   = omp_get_wtime()

      tree_set%bsize = bsize
      tree_set%dtype = d_type
      tree_set%vtype = v_type
      tree_set%ndim = n_dim
      tree_set%n_thread = n_thread
      !tree_set%np_dmax_tag = -1
      !tree_set%np_mass_tag = 1
      root = js_kdtree_mktree(pos, mm, orgind, tree_set)!bsize, d_type, v_type, n_dim, n_thread)

      time(2)   = omp_get_wtime()
      !!-----
      !! GET LEAF ONLY
      !!-----
      n_leaf    = js_kdtree_getleafnum(root)
      ALLOCATE(lf(1:n_leaf))

      CALL js_kdtree_getleaf(root, lf)

      !!-----
      !! COMPUTE POTENTIAL
      !!-----

      !!----- FOR DS TEST
      !l = 0
      !!$OMP PARALLEL DO default(shared) &
      !!$OMP & private(j, dx, k) schedule(static) reduction(+:l)
      !DO i=1, n_ptcl
      !  IF(pot(i) .GT. 0) CYCLE
      !  k       = omp_get_thread_num()
      !  IF(k .EQ. 1) PRINT *, l
      !  DO j=1, n_ptcl
      !    IF(i .EQ. j)CYCLE
      !    dx = 0.
      !    DO k=1, n_dim
      !      dx = dx + ( pos(i,k) - pos(j,k) )**2
      !    ENDDO
      !    dx = dx**0.5

      !    pot(i) = pot(i) + (-Gconst) * mm(j) / dx
      !  ENDDO
      !  l = l + 1
      !ENDDO
      !!$OMP END PARALLEL DO
      !CALL QUICKSORT(orgind, recind, SIZE(orgind), 1, SIZE(orgind))

      !DO i=1, n_dim
      !  pos(:,i) = pos(recind,i)
      !ENDDO
      !mm = mm(recind)
      !pot = pot(recind)
      !RETURN

      IF(e_type .EQ. 1) THEN 
      !!----- Considering Particle position
        !$OMP PARALLEL DO default(shared) &
        !$OMP & private(bs, be, dx, j, k, l) &
        !$OMP & schedule(static)
        DO i=1, n_ptcl
          pot(i) = 0.
          force(i) = 0.

          IF(p_type .EQ. 0) THEN !! mesh to mesh only

            DO j=1, n_leaf
              bs = lf(j)%bstart
              be = lf(j)%bend

              dx = 0.
              DO k=1, n_dim
                dx = dx + ( lf(j)%cen(k) - pos(i,k) )**2
              ENDDO
              dx = dx**0.5

              IF(i .GE. bs .AND. i .LE. be) THEN
                pot(i) = pot(i) + (-Gconst) * (lf(j)%mass - mm(i))/dx
                force(i) = force(i) + (-Gconst) * (lf(j)%mass - mm(i))/dx**2
              ELSE
                pot(i) = pot(i) + (-Gconst) * lf(j)%mass/dx
                force(i) = force(i) + (-Gconst) * lf(j)%mass/dx**2
              ENDIF
            ENDDO

          ENDIF

          IF(p_type .EQ. 1) THEN !! particle to mesh
            DO j=1, n_leaf
              bs = lf(j)%bstart
              be = lf(j)%bend

              IF(i .GE. bs .AND. i .LE. be) THEN
                DO k=bs, be
                  IF(k .EQ. i) CYCLE

                  dx = 0.
                  DO l=1, n_dim
                    dx = dx + ( pos(i,l) - pos(k,l) )**2
                  ENDDO
                  dx = dx**0.5

                  pot(i) = pot(i) + (-Gconst) * mm(k)/dx
                  force(i) = force(i) + (-Gconst) * mm(k)/dx**2
                ENDDO
              ELSE
                dx = 0.
                DO k=1, n_dim
                  dx = dx + ( pos(i,k) - lf(j)%cen(k) )**2
                ENDDO
                dx = dx**0.5

                pot(i) = pot(i) + (-Gconst) * lf(j)%mass/dx
                force(i) = force(i) + (-Gconst) * lf(j)%mass/dx**2
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF

      !IF(e_type .EQ. 0) THEN
      !  !----- Particle position is fixed to the center of containing
      !  !        node - cost efficient
      !  !j, dummy, dx, bs, be, k
      !  !$OMP PARALLEL DO default(shared) &
      !  !$OMP & private(dummy, j, k, l, bs, be, dx) &
      !  !$OMP & schedule(static)
      !  DO i=1, n_leaf

      !    ! Monopole term
      !    dummy = 0.
      !    DO j=1, n_leaf
      !      IF(i .EQ. j) CYCLE
      !      dx = js_getpt_ft_d3d(lf, i, j, n_dim)
      !      IF(dx .EQ. 0) CYCLE
      !      dummy = dummy + (-Gconst) * lf(j)%mass / dx
      !    ENDDO

      !    ! Direct term
      !    bs      = lf(i)%bstart
      !    be      = lf(i)%bend

      !    DO j=bs, be
      !      pot(j) = pot(j) + dummy

      !      IF(p_type .EQ. 0) THEN !! mesh to mesh only
      !        dx = 0.
      !        DO k=1, n_dim
      !          dx = dx + ( lf(i)%cen(k) - pos(j,k) )**2
      !        ENDDO
      !        dx = dx**0.5
      !        IF(dx .EQ. 0) CYCLE
      !        pot(j) = pot(j) + (-Gconst) * (lf(i)%mass - mm(j)) / dx

      !      ENDIF

      !      IF(p_type .EQ. 1) THEN !! particle to mesh
      !        DO k=bs, be
      !          IF(j .EQ. k) CYCLE

      !          dx = 0.
      !          DO l=1, n_dim
      !            dx = dx + ( pos(j,l) - pos(k,l) )**2
      !          ENDDO
      !          dx = dx**0.5
      !          IF(dx .EQ. 0) CYCLE
      !          pot(j) = pot(j) + (-Gconst) * mm(k) / dx
      !        ENDDO
      !      ENDIF
      !    ENDDO
      !  ENDDO
      !  !$OMP END PARALLEL DO
      !ENDIF

      IF(e_type .EQ. 0) THEN
        !----- Particle position is fixed to the center of containing
        !        node - cost efficient
        !j, dummy, dx, bs, be, k
        !$OMP PARALLEL DO default(shared) &
        !$OMP & private(dummy, dummy_v, dummy_vf, bnd, j, k, l, bs, be, dx) &
        !$OMP & schedule(static)
        DO i=1, n_leaf

          ! Monopole term
          dummy_v = 0.
          dummy_vf = 0.
          bnd   = js_getpt_ft_getbnd(lf(i))

          DO j=1, n_leaf
            IF(i .EQ. j) CYCLE

            DO k=1, 8
              dx = js_getpt_ft_dist(lf(j)%cen(1:3),bnd(k,:))
              IF(dx .EQ. 0) CYCLE
              dummy_v(k) = dummy_v(k) + (-Gconst) * lf(j)%mass / dx
              dummy_vf(k) = dummy_vf(k) + (-Gconst) * lf(j)%mass / dx**2
            ENDDO
          ENDDO

          ! Direct term
          bs      = lf(i)%bstart
          be      = lf(i)%bend

          DO j=bs, be
            dummy  = js_getpt_ft_interpole(bnd, pos(j,1:3), dummy_v)
            pot(j) = pot(j) + dummy


            dummy  = js_getpt_ft_interpole(bnd, pos(j,1:3), dummy_vf)
            force(j) = force(j) + dummy

            IF(p_type .EQ. 0) THEN !! mesh to mesh only
              dx = 0.
              DO k=1, n_dim
                dx = dx + ( lf(i)%cen(k) - pos(j,k) )**2
              ENDDO
              dx = dx**0.5
              IF(dx .EQ. 0) CYCLE
              pot(j) = pot(j) + (-Gconst) * (lf(i)%mass - mm(j)) / dx
              force(j) = force(j) + (-Gconst) * (lf(i)%mass - mm(j)) / dx**2
            ENDIF

            IF(p_type .EQ. 1) THEN !! particle to mesh
              DO k=bs, be
                IF(j .EQ. k) CYCLE

                dx = 0.
                DO l=1, n_dim
                  dx = dx + ( pos(j,l) - pos(k,l) )**2
                ENDDO
                dx = dx**0.5
                IF(dx .EQ. 0) CYCLE
                pot(j) = pot(j) + (-Gconst) * mm(k) / dx
                force(j) = force(j) + (-Gconst) * mm(k) / dx**2
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF

      time(3)   = omp_get_wtime()

      !!-----
      !! REORDERING
      !!-----

      CALL QUICKSORT(orgind, recind, SIZE(orgind), 1, SIZE(orgind))

      !DO i=1, n_dim
      !  pos(:,i) = pos(recind,i)
      !ENDDO
      !mm = mm(recind)
      pot = pot(recind)
      force = force(recind)

      PRINT *, '%123123---------------'
      PRINT *, '        Wall-clock time Report'
      PRINT *, '        Tree in ', time(2) - time(1)
      PRINT *, '        Pot in ', time(3) - time(2)
      PRINT *, '%123123---------------'
      DEALLOCATE(lf)
      DEALLOCATE(root)
      !DEALLOCATE(part)
      RETURN
CONTAINS

      FUNCTION js_getpt_ft_d3d(node, i, j, ndim) RESULT(dx)
        USE js_kdtree
        IMPLICIT NONE
        TYPE(nodetype), DIMENSION(:) :: node
        INTEGER(KIND=4) i, j, k, ndim
        REAL(KIND=8) dx

        dx = 0.
        DO k=1, ndim
          dx = dx + (node(i)%cen(k) - node(j)%cen(k))**2
        ENDDO
        dx = dx**0.5
        RETURN
      END FUNCTION js_getpt_ft_d3d

      FUNCTION js_getpt_ft_dist(x0, x1) RESULT(dd)
        IMPLICIT NONE
        REAL(KIND=8) x0(3), x1(3)
        REAL(KIND=8) dd
        INTEGER(KIND=4) i

        dd = 0.
        DO i=1, 3
          dd = dd + ( x0(i) - x1(i) )**2
        ENDDO
        dd = dd**0.5

        RETURN
      END FUNCTION js_getpt_ft_dist

      FUNCTION js_getpt_ft_getbnd(node) RESULT(bnd)
        USE js_kdtree
        IMPLICIT NONE
        TYPE(nodetype) node
        REAL(KIND=8) bnd(8,3)
        INTEGER(KIND=4) i, j

        DO i=1, 2
          j = 4*(i-1)
          bnd(j+1,1) = node%bnd(1,1)
          bnd(j+1,2) = node%bnd(2,1)
          bnd(j+1,3) = node%bnd(3,i)

          bnd(j+2,1) = node%bnd(1,1)
          bnd(j+2,2) = node%bnd(2,2)
          bnd(j+2,3) = node%bnd(3,i)

          bnd(j+3,1) = node%bnd(1,2)
          bnd(j+3,2) = node%bnd(2,2)
          bnd(j+3,3) = node%bnd(3,i)

          bnd(j+4,1) = node%bnd(1,2)
          bnd(j+4,2) = node%bnd(2,1)
          bnd(j+4,3) = node%bnd(3,i)
        ENDDO
        RETURN
      END FUNCTION js_getpt_ft_getbnd

      FUNCTION js_getpt_ft_interpole(bnd, pos, dummy_v) RESULT(p)
        IMPLICIT NONE
        REAL(KIND=8) bnd(8,3), pos(3), dummy_v(8), p
        REAL(KIND=8) v1, v2, v0
        REAL(KIND=8) w1, w2, w0

        v1      = js_getpt_ft_lint(pos(2), bnd(1,2), bnd(2,2), dummy_v(1), dummy_v(2))
        v2      = js_getpt_ft_lint(pos(2), bnd(4,2), bnd(3,2), dummy_v(4), dummy_v(3))
        v0      = js_getpt_ft_lint(pos(1), bnd(1,1), bnd(3,1), v1, v2)

        w1      = js_getpt_ft_lint(pos(2), bnd(5,2), bnd(6,2), dummy_v(5), dummy_v(6))
        w2      = js_getpt_ft_lint(pos(2), bnd(8,2), bnd(7,2), dummy_v(8), dummy_v(7))
        w0      = js_getpt_ft_lint(pos(1), bnd(5,1), bnd(7,1), w1, w2)

        p       = js_getpt_ft_lint(pos(3), bnd(1,3), bnd(5,3), v0, w0)
        RETURN
      END FUNCTION js_getpt_ft_interpole
        
      FUNCTION js_getpt_ft_lint(x, x0, x1, v0, v1) RESULT(val)
        IMPLICIT NONE
        REAL(KIND=8) x, x0, x1, v0, v1, val
        REAL(KIND=8) dx

        dx = (x-x0) / (x1-x0)
        dx = MAX(dx,0d0)
        val= dx*v1 + (1d0 - dx)*v0
        RETURN
      END FUNCTION js_getpt_ft_lint


!! QUICK SORT
!!!!!
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
      recursive subroutine quicksort(a, b, nn, first, last)
        implicit none
        integer*4  a(nn), b(nn), x, t, n
        integer first, last
        integer i, j, nn

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
           do while (a(i) < x)
              i=i+1
           end do
           do while (x < a(j))
              j=j-1
           end do
           IF(i .GE. j) exit
           t = a(i);  a(i) = a(j);  a(j) = t
           n = b(i); b(i) = b(j); b(j) = n
           i=i+1
           j=j-1
        end do
        if (first < i-1) call quicksort(a, b, nn, first, i-1)
        if (j+1 < last)  call quicksort(a, b, nn, j+1, last)
        return
      end subroutine quicksort
      END SUBROUTINE js_getpt_ft
