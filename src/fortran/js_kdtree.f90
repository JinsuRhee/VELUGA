MODULE js_kdtree
!234567

      TYPE nodetype
        INTEGER(KIND=4) :: id, bstart, bend, ncount
        INTEGER(KIND=4) :: splitdim, splitind, leaf
        REAL(KIND=8) :: bnd(6,2), cen(6), mass, splitval
        INTEGER(KIND=4) :: numnode, level
        INTEGER(KIND=4) :: parent, sibling, left, right
        REAL(KIND=8) :: dmax
      END TYPE nodetype

      TYPE infotype
        INTEGER(KIND=4) :: bsize, ndim, npart
        INTEGER(KIND=4) :: dtype, vtype
        INTEGER(KIND=4) :: omp_tag, n_thread, dmax_tag
      END TYPE infotype
CONTAINS

!!--------------------------------------------------
!! Main Routine
!!--------------------------------------------------
!234567
      FUNCTION js_kdtree_mktree(pos, mm, orgind, &
        info) RESULT(root)

        USE omp_lib
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:) :: pos
        REAL(KIND=8), DIMENSION(:) :: mm
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: posdum
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: mmdum
        INTEGER(KIND=4), DIMENSION(:) :: orgind
        !INTEGER(KIND=4), INTENT(IN) :: bsize, dtype, vtype, ndim, n_thread, n_thread2
        INTEGER(KIND=4) :: act_thread, bsize_p
        INTEGER(KIND=4) :: npart, i, k
        INTEGER(KIND=4) :: numnode, level, bstart, bend
        TYPE(nodetype), DIMENSION(:), ALLOCATABLE :: root, root_dum, root_p, root_pl, root_t
        TYPE(nodetype), DIMENSION(:,:), ALLOCATABLE :: root_2d
        !TYPE(dat), DIMENSION(:), INTENT(inout) :: part
        TYPE(infotype) info, info2
        INTEGER(KIND=4) :: bs, be, bs0, be0, n_ini, n_aft
        !TYPE(dat), DIMENSION(:), ALLOCATABLE :: partdum
        INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: inddum, nums
        INTEGER(KIND=4) idoffset, ind0, ind1, lind, rind, indoffset
        INTEGER(KIND=4) memdebug, mem_nn

        !! Initialize
        npart   = SIZE(mm)
        !info%bsize      = bsize
        !info%ndim       = ndim
        !info%npart      = npart
        !info%dtype      = dtype
        !info%vtype      = vtype
        !info%omp_tag    = -1
        !info%n_thread   = n_thread
        !info%dmax_tag   = -1
        !! Set Threads
        k = 0.
        DO
          IF(2.**k .LT. info%n_thread+0.1) THEN
            k = k+1
          ELSE
            act_thread = 2.**(k-1)
            EXIT
          ENDIF
        ENDDO

        !! Make initial nodes for parallelization
        IF(ALLOCATED(root)) DEALLOCATE(root)
        ALLOCATE(root(1:npart))

        info2 = info
        info2%bsize = npart / act_thread + 1
        info2%dtype = 0
        info2%vtype = 0
        info2%omp_tag = 1

        DO i=1, info%ndim
          root(1)%bnd(i,1) = MINVAL(pos(:,i))
          root(1)%bnd(i,2) = MAXVAL(pos(:,i))
        ENDDO

        numnode = 0
        level = 1
        bstart = 1
        bend = npart
        CALL js_kdtree_buildnode(root, pos, mm, orgind, info2, &
                numnode, bstart, bend, level)

        IF(info2%bsize .LE. info%bsize) THEN
          !!----- No need to make further son nodes
          RETURN
        ENDIF

        IF(ALLOCATED(root_p)) DEALLOCATE(root_p)
        ALLOCATE(root_p(1:numnode))
        root_p = root(1:numnode)
        DEALLOCATE(root)

        !! Leaf for seed
        IF(ALLOCATED(root_pl)) DEALLOCATE(root_pl)
        ALLOCATE(root_pl(1:act_thread))


        k = 1
        DO i=1, numnode
          IF(root_p(i)%leaf .EQ. 1) THEN
            root_pl(k) = root_p(i)
            k = k+1
          ENDIF
        ENDDO
        n_ini      = numnode

        IF(k-1 .NE. act_thread) STOP
        !!----- Build Node
        ! bs, be, partdum, inddum, numnode, level
        ! bs0, be0
        ! n_aft
        IF (ALLOCATED(root_2d)) DEALLOCATE(root_2d)
        ALLOCATE(root_2d(1:npart/act_thread, 1:act_thread))
        
        IF (ALLOCATED(nums)) DEALLOCATE(nums)
        ALLOCATE(nums(1:act_thread), stat=memdebug)

        nums = 0
        n_aft = 0

        !$OMP PARALLEL DO default(shared) &
        !$OMP & private(bs, be, posdum, mmdum, inddum, numnode) &
        !$OMP & private(level, bs0, be0) &
        !$OMP & reduction(+:n_aft)
        DO i=1, act_thread
          bs = root_pl(i)%bstart
          be = root_pl(i)%bend

          ALLOCATE(posdum(1:be-bs+1,1:info%ndim))
          ALLOCATE(mmdum(1:be-bs+1))
          ALLOCATE(inddum(1:be-bs+1))
          posdum = pos(bs:be,:)
          mmdum   = mm(bs:be)
          inddum = orgind(bs:be)


          !DO k=1, ndim
          !  root_2d(1,i)%bnd(k,1) = js_kdtree_min(partdum%pos(k))
          !  root_2d(1,i)%bnd(k,2) = js_kdtree_max(partdum%pos(k))
          !ENDDO
          root_2d(1,i)%bnd = root_pl(i)%bnd

          numnode = 0
          level = 1

          bs0 = 1
          be0 = be-bs+1
          CALL js_kdtree_buildnode(root_2d(:,i), posdum, mmdum, inddum, info, &
            numnode, bs0, be0, level)
          nums(i) = numnode
          n_aft = n_aft + nums(i) - 1
          pos(bs:be,:) = posdum
          mm(bs:be) = mmdum
          orgind(bs:be) = inddum
          DEALLOCATE(posdum, mmdum, inddum)
        ENDDO
        !$OMP END PARALLEL DO

        !!--Merge
        ALLOCATE(root(1:n_ini+n_aft))


        !!!!----- Initial node
        root(1:n_ini) = root_p
        idoffset = root(n_ini)%id - 1
        indoffset = 0
        ind0    = n_ini+1

        DO i=1, act_thread
          ind1 = ind0 + nums(i) - 2

          lind  = root_2d(1,i)%left
          rind  = root_2d(1,i)%right

          root_2d(1:nums(i),i)%id = root_2d(1:nums(i),i)%id + idoffset
          root_2d(1:nums(i),i)%numnode = root_2d(1:nums(i),i)%numnode + idoffset
          root_2d(1:nums(i),i)%parent = root_2d(1:nums(i),i)%parent + idoffset
          root_2d(1:nums(i),i)%left = root_2d(1:nums(i),i)%left + idoffset
          root_2d(1:nums(i),i)%right = root_2d(1:nums(i),i)%right + idoffset
          root_2d(1:nums(i),i)%sibling = root_2d(1:nums(i),i)%sibling + idoffset
          root_2d(1:nums(i),i)%bstart = root_2d(1:nums(i),i)%bstart + indoffset
          root_2d(1:nums(i),i)%bend = root_2d(1:nums(i),i)%bend + indoffset

          root_2d(lind,i)%parent = root_pl(i)%id
          root_2d(rind,i)%parent = root_pl(i)%id

          root(root_pl(i)%id)%left = root_2d(lind,i)%id
          root(root_pl(i)%id)%right = root_2d(rind,i)%id
          root(root_pl(i)%id)%leaf = -1

          root(ind0:ind1) = root_2d(2:nums(i),i)

          ind0  = ind1 + 1
          idoffset = idoffset + nums(i)-1
          indoffset = root_2d(nums(i),i)%bend

        ENDDO


        DEALLOCATE(root_p, root_pl, root_2d, nums)
        root(1)%numnode = n_ini + n_aft

      !  PRINT *, ind1, n_ini, n_aft, n_ini+n_aft
      !  PRINT *, nums
      !  PRINT *, npart
      !  DO i=1, n_ini+n_aft
      !    PRINT *, 'ID = ', root(i)%id,'     i = ', i
      !    PRINT *, 'IND = ', root(i)%bstart, ' - ', root(i)%bend
      !    PRINT *, '    Par = ', root(i)%parent, 'L/R = ', root(i)%left, root(i)%right
      !    PRINT *,'     Sib = ', root(i)%sibling
      !    PRINT *,'     Leaf = ', root(i)%leaf
      !    PRINT *,'  ' 
      !  ENDDO
      !  DEALLOCATE(root_p, root_pl, root_2d, nums)
      !  STOP

!        ALLOCATE(root(1:npart))
!        ALLOCATE(root_dum(1:npart))
!
!        !!----- ROOT NODE INFORMATION
!        DO i=1, ndim
!          root(1)%bnd(i,1) = js_kdtree_min(part%pos(i))
!          root(1)%bnd(i,2) = js_kdtree_max(part%pos(i))
!        ENDDO
!
!
!        !!----- Build NODE
!        numnode = 0
!        level   = 1
!        bstart  = 1
!        bend    = npart
!
!        !!$omp parallel default(shared) private(bstart, bend, level)
!        CALL js_kdtree_buildnode(root, part, orgind, info, &
!                numnode, bstart, bend, level)
!
!        !!$omp end parallel
!        !!----- NEW INDEX
!        !DO i=1, ndim
!        !  pos(:,i) = part%pos(i)
!        !ENDDO
!        !mm(:)   = part%mm
!
!        !DEALLOCATE(part)
!
!        !!----- DEALLOCATE
!        root_dum        = root
!        DEALLOCATE(root)
!        ALLOCATE(root(1:numnode))
!        root    = root_dum(1:numnode)
!        DEALLOCATE(root_dum)
!        root(1)%numnode = numnode
      END FUNCTION js_kdtree_mktree

!!--------------------------------------------------
!! BUILD NODE
!!--------------------------------------------------
      RECURSIVE SUBROUTINE js_kdtree_buildnode(node, pos, mm, orgind, info, &
                      numnode, bstart, bend, level)
        USE omp_lib
        IMPLICIT NONE
        
        TYPE(nodetype), DIMENSION(:) :: node
        !TYPE(dat), DIMENSION(:) :: part
        REAL(KIND=8), DIMENSION(:,:) :: pos
        REAL(KIND=8), DIMENSION(:) :: mm
        !TYPE(dat), DIMENSION(:), ALLOCATABLE :: partdum
        TYPE(infotype), INTENT(IN) :: info
        INTEGER(KIND=4), DIMENSION(:) :: orgind
        INTEGER(KIND=4) :: numnode, bstart, bend, level, level_up, nid, bs, be

        INTEGER(KIND=4) i, j, k
        !REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: posdum
        !REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mdum
        !REAL(KIND=8) debug_time(10)     !123123

        !debug_time(1) = omp_get_wtime()    !123123
        !!----- SET NODE PROPs
        numnode = numnode + 1
        node(numnode)%id      = numnode
        node(numnode)%bstart  = bstart
        node(numnode)%bend    = bend
        node(numnode)%ncount  = bend-bstart+1
        node(numnode)%level   = level
        node(numnode)%leaf    = -1

        !debug_time(1) = omp_get_wtime()    !123123
        !!----- FIRST EXTRACT TARGET

        DO i=1, info%ndim
          IF(SUM(mm(bstart:bend)) .GT. 0) THEN
            node(numnode)%cen(i) = SUM(pos(bstart:bend,i) * mm(bstart:bend)) / SUM(mm(bstart:bend))
          ELSE
            node(numnode)%cen(i) = SUM(pos(bstart:bend,i)) / (bend-bstart + 1.)
          ENDIF
        ENDDO

        ! for dmax
        IF(info%dmax_tag .EQ. 1) node(numnode)%dmax = js_kdtree_dmax(pos(bstart:bend,:), node(numnode)%cen, info%ndim)
        !!----- IS LEAF (some props are only calculated for leafs)
        IF(node(numnode)%ncount .LE. info%bsize) THEN
          !for mass
          node(numnode)%mass  = SUM(mm(bstart:bend))!js_kdtree_total(info, mdum)
          node(numnode)%leaf = 1
          RETURN
        ENDIF

        !!----- DETERMINE SPLITDIM
        node(numnode)%splitdim = js_kdtree_sdim(pos(bstart:bend,:), info)
                
        !!----- DETERMIN SPLITVALUE
        CALL js_kdtree_sval(node(numnode), pos, mm, bstart, bend, orgind, info)
        
        !debug_time(5) = omp_get_wtime()    !123123
        !PRINT *, '      value made in :', debug_time(5)-debug_time(4)

        
        !!----- BUILD Sons

        nid     = numnode

        level_up = level + 1
        node(numnode)%left      = numnode + 1
        node(numnode+1)%parent  = nid
        node(node(nid)%left)%bnd        = node(nid)%bnd
        node(node(nid)%left)%bnd( (node(nid)%splitdim), 2) = node(nid)%splitval

        bs  = node(nid)%bstart
        be  = node(nid)%splitind - 1
        
        CALL js_kdtree_buildnode(node, pos, mm, orgind, info, &
                numnode, bs, be, level_up)

        level_up = level+1
        node(nid)%right = numnode + 1
        node(numnode+1)%parent = nid
        node(node(nid)%right)%bnd       = node(nid)%bnd
        node(node(nid)%right)%bnd( (node(nid)%splitdim), 1) = node(nid)%splitval

        bs    = node(nid)%splitind
        be    = node(nid)%bend
        
        CALL js_kdtree_buildnode(node, pos, mm, orgind, info, &
                numnode, bs, be, level_up)

        node(node(nid)%left)%sibling = node(nid)%right
        node(node(nid)%right)%sibling = node(nid)%left

        !!----- Set SONs
        !node%left       => left
        !node%right      => right
        !node%left       => lnode
        !node%right      => rnode

        !left%parent     => node
        !right%parent    => node

        !left%sibling    => right
        !right%sibling   => left

        !lnode%bnd     = node%bnd
        !rnode%bnd    = node%bnd

        !lnode%bnd(node%splitdim, 2)        = node%splitval
        !rnode%bnd(node%splitdim, 1)       = node%splitval

        !node(node(nid)%left)%bnd        = node(nid)%bnd
        !node(node(nid)%right)%bnd       = node(nid)%bnd

        !node(node(nid)%left)%bnd( (node(nid)%splitdim), 2) = node(nid)%splitval
        !node(node(nid)%right)%bnd( (node(nid)%splitdim), 1) = node(nid)%splitval

        RETURN
      END SUBROUTINE js_kdtree_buildnode

!!--------------------------------------------------
!! TREE SHAPE FTNs
!!--------------------------------------------------
      !! SPLIT VALUE
      SUBROUTINE js_kdtree_sval(node, pos, mass, bstart, bend, orgind, info)
        IMPLICIT NONE
        !TYPE(dat), DIMENSION(:) :: part
        INTEGER(KIND=4) bstart, bend
        INTEGER(KIND=4), DIMENSION(:) :: orgind
        TYPE(infotype) info
        TYPE(nodetype) node
        REAL(KIND=8), DIMENSION(:,:) :: pos
        REAL(KIND=8), DIMENSION(:) :: mass
        !REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: :
        !REAL(KIND=8), DIMENSION(:), INTENT(IN) :: massdum

        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tmp
        INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ind
        INTEGER(KIND=4) i, j, k, sind

        IF(info%vtype .EQ. 0) THEN ! BALANCED
          ALLOCATE(tmp(1:(bend-bstart+1)))
          ALLOCATE(ind(1:(bend-bstart+1)))

          tmp     = pos(bstart:bend,node%splitdim)
          DO i=1, SIZE(tmp)
            ind(i) = i
          ENDDO

          CALL js_kdtree_mediansort(tmp, ind, sind)
          ind = ind + bstart - 1
          sind = sind + bstart - 1
          pos(bstart:bend,:) = pos(ind,:)
          mass(bstart:bend) = mass(ind)
          orgind(bstart:bend) = orgind(ind) 

          node%splitval = pos(sind,node%splitdim)
          node%splitind = sind

          DEALLOCATE(tmp)
          DEALLOCATE(ind)
        ELSE
          PRINT *, 'not implemented yet'
          STOP
        ENDIF
      END SUBROUTINE js_kdtree_sval

      !! SPLIT DIMENSION
      INTEGER(KIND=4) FUNCTION js_kdtree_sdim(pos, info)
        IMPLICIT NONE
        !TYPE(dat), DIMENSION(:), INTENT(IN) :: part
        REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: pos
        TYPE(infotype) :: info
        !INTEGER(KIND=4) :: bstart, bend
        INTEGER(KIND=4) :: i, j, k
        REAL(KIND=8) :: dx, dx2

        js_kdtree_sdim = 1
        IF(info%dtype .EQ. 0) THEN ! MAX RANGE
          dx = 0.
          DO i=1,info%ndim
            !dx2 = js_kdtree_max(part%pos(i)) - js_kdtree_min(part%pos(i))
            dx2 = MAXVAL(pos(:,i)) - MINVAL(pos(:,i))
            IF(dx2 .GE. dx) THEN
              dx = dx2
              js_kdtree_sdim = i
            ENDIF
          ENDDO
        ELSE 
          PRINT *, 'not implemented yet'
          STOP
        ENDIF
      END FUNCTION js_kdtree_sdim

        
!!--------------------------------------------------
!! SIMPLE FTNs
!!--------------------------------------------------
      ! MIN
      REAL(KIND=8) FUNCTION js_kdtree_min(xx)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT(IN) ::xx
        INTEGER(KIND=4) i

        js_kdtree_min = xx(1)
        DO i=2, SIZE(xx)
          IF(xx(i) .LE. js_kdtree_min) js_kdtree_min = xx(i)
        ENDDO
      END FUNCTION js_kdtree_min

      !MAX
      REAL(KIND=8) FUNCTION js_kdtree_max(xx)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xx
        INTEGER(KIND=4) i

        js_kdtree_max = xx(1)
        DO i=2, SIZE(xx)
          IF(xx(i) .GE. js_kdtree_max) js_kdtree_max = xx(i)
        ENDDO
      END FUNCTION js_kdtree_max

      !TOTAL
      REAL(KIND=8) FUNCTION js_kdtree_total(info, xx)
        USE omp_lib
        IMPLICIT NONE
        TYPE(infotype), INTENT(IN) :: info
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xx
        INTEGER(KIND=4) i

        js_kdtree_total = 0.
        IF(info%omp_tag .EQ. 1) THEN
          CALL omp_set_num_threads(info%n_thread)
          !$OMP PARALLEL DO default(shared) &
          !$OMP & reduction(+:js_kdtree_total)
          DO i=1, SIZE(xx)
            js_kdtree_total = js_kdtree_total + xx(i)
            IF(i .EQ. 1) PRINT *, omp_get_num_threads()
          ENDDO
          !$OMP END PARALLEL DO
          PRINT *, js_kdtree_total
        ELSE
          DO i=1, SIZE(xx)
            js_kdtree_total = js_kdtree_total + xx(i)
          ENDDO
        ENDIF
      END FUNCTION js_kdtree_total

      !Compute metric
      REAL(KIND=8) FUNCTION js_kdtree_dmax (pos, yy, ndim)
        IMPLICIT NONE
        !TYPE(dat), DIMENSION(:), INTENT(IN) :: part
        REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: pos
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: yy
        INTEGER(KIND=4) ndim, npart, i
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ddum

        npart = SIZE(pos(:,1))
        ALLOCATE(ddum(1:npart))

        ddum = 0.

        DO i=1, ndim
          ddum = ddum + (pos(:,i) - yy(i))**2
        ENDDO
        ddum = SQRT(ddum)

        js_kdtree_dmax = MAXVAL(ddum)!js_kdtree_max(ddum)
        DEALLOCATE(ddum)
      END FUNCTION js_kdtree_dmax

      !SWAP FTN
      SUBROUTINE js_kdtree_swap_l(arr, i, j)
        IMPLICIT NONE
        INTEGER(KIND=4) i, j
        INTEGER(KIND=4), DIMENSION(:) :: arr
        INTEGER(KIND=4) dum

        dum = arr(i)
        arr(i) = arr(j)
        arr(j) = dum
      END SUBROUTINE js_kdtree_swap_l
      SUBROUTINE js_kdtree_swap_d(arr, i, j)
        IMPLICIT NONE
        INTEGER(KIND=4) i, j
        REAL(KIND=8), DIMENSION(:) :: arr
        REAL(KIND=8) dum

        dum = arr(i)
        arr(i) = arr(j)
        arr(j) = dum
      END SUBROUTINE js_kdtree_swap_d

      !MEDIAN SORT
      SUBROUTINE js_kdtree_mediansort(tmp, ind, sind)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:) :: tmp
        INTEGER(KIND=4), DIMENSION(:) :: ind

        INTEGER(KIND=4) nn, i0, i1, j0, j1, k, pivot_ind, sind
        REAL(KIND=8) x, pivot

        i0 = 1
        j0 = SIZE(tmp)
        k  = i0 + (j0-i0 + 1)/2

        DO WHILE(i0 .LT. j0)
          x = tmp(k)
          CALL js_kdtree_swap_d(tmp, j0, k)
          CALL js_kdtree_swap_l(ind, j0, k)

          pivot = tmp(k)
          i1    = i0 - 1
          j1    = j0

          DO WHILE (.TRUE.)
            DO WHILE(i1 .LT. j1)
              i1 = i1 + 1
              IF(tmp(i1) .GE. x) EXIT
            ENDDO

            DO WHILE(i1 .LT. j1)
              j1 = j1 - 1
              IF(tmp(j1) .LE. x) EXIT
            ENDDO

            CALL js_kdtree_swap_d(tmp, i1, j1)
            CALL js_kdtree_swap_l(ind, i1, j1)
            pivot = tmp(j1)
            pivot_ind = ind(j1)
            IF( j1 .LE. i1) EXIT
          ENDDO

          tmp(j1)     = tmp(i1)
          ind(j1)       = ind(i1)

          tmp(i1)     = tmp(j0)
          ind(j1)       = ind(j0)

          tmp(j0)     = pivot
          ind(j0)       = pivot_ind

          IF(i1 .GE. k) j0 = i1 - 1
          IF(i1 .LE. k) i0 = i1 + 1
        ENDDO
        sind = k
      END SUBROUTINE js_kdtree_mediansort

      ! get number of leaf nodes
      RECURSIVE SUBROUTINE js_kdtree_getleafnum_wktree(node, nn, nid)
        TYPE(nodetype), DIMENSION(:) :: node
        INTEGER(KIND=4) nn, left, right, nid

        IF(node(nid)%leaf .GT. 0) THEN
          nn = nn + 1
          RETURN
        ELSE
          left = node(nid)%left
          CALL js_kdtree_getleafnum_wktree(node, nn, left)

          right = node(nid)%right
          CALL js_kdtree_getleafnum_wktree(node, nn, right)
        ENDIF
      END SUBROUTINE js_kdtree_getleafnum_wktree
      FUNCTION js_kdtree_getleafnum(node) RESULT(n_leaf)
        TYPE(nodetype), DIMENSION(:) :: node
        INTEGER(KIND=4) n_leaf, nid
     
        n_leaf = 0
        nid = 1 
        CALL js_kdtree_getleafnum_wktree(node, n_leaf, nid)

        RETURN
      END FUNCTION js_kdtree_getleafnum

      ! GET LEAF ONLY
      RECURSIVE SUBROUTINE js_kdtree_getleaf_fill(node, lf, nid, lid)
        TYPE(nodetype), DIMENSION(:) :: node, lf
        INTEGER(KIND=4) nid, lid, left, right

        IF(node(nid)%leaf .GT. 0) THEN
          lf(lid) = node(nid)
          lid = lid + 1
        ELSE
          left = node(nid)%left
          right = node(nid)%right

          CALL js_kdtree_getleaf_fill(node, lf, left, lid)
          CALL js_kdtree_getleaf_fill(node, lf, right, lid)
        ENDIF
      END SUBROUTINE js_kdtree_getleaf_fill
      SUBROUTINE js_kdtree_getleaf(root, lf)
        TYPE(nodetype), DIMENSION(:) :: root, lf
        INTEGER(KIND=4) nid, lid

        nid = 1
        lid = 1
        CALL js_kdtree_getleaf_fill(root, lf, nid, lid)

      END SUBROUTINE js_kdtree_getleaf

      FUNCTION js_kdtree_nodeskip(node, pos, r, ndim) RESULT(skip)
        IMPLICIT NONE
        LOGICAL skip
        TYPE(nodetype), INTENT(IN) :: node
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: pos
        REAL(KIND=8), INTENT(IN) :: r
        INTEGER(KIND=4), INTENT(IN) :: ndim

        REAL(KIND=8) d_max, d_center, Dmin, Dmax
        INTEGER(KIND=4) i
        skip = .false.

        d_max = node%dmax
        d_center = 0.
        DO i=1, ndim
          d_center = d_center + (node%cen(i) - pos(i))**2
        ENDDO
        d_center = SQRT(d_center)


        IF(d_max .GE. d_center) THEN
          skip = .false.
          RETURN
        ELSE
          Dmin = d_center - d_max
          IF(Dmin .GT. r) THEN
            skip = .true.
            RETURN
          ELSE
            skip = .false.
            RETURN
          ENDIF
        ENDIF
      END FUNCTION js_kdtree_nodeskip


END MODULE
