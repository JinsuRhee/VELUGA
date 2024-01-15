!1234567        
       SUBROUTINE get_contam(larr, darr, dir_raw, xc, yc, zc, aperture, &
                       xp, mp, conf_n, conf_m)

       USE omp_lib
       USE js_kdtree

!!-----
!! GLOBAL VARIABLES
!-----
       IMPLICIT NONE

       REAL(KIND=8) darr(20)
       INTEGER(KIND=4) larr(20)

       REAL(KIND=8) xc(larr(1)), yc(larr(1)), zc(larr(1))
       REAL(KIND=8) aperture(larr(1),larr(3))

       REAL(KIND=8) xp(larr(2), 3), mp(larr(2))
       REAL(KIND=8) conf_n(larr(1), larr(3))
       REAL(KIND=8) conf_m(larr(1), larr(3))
       !REAL(KIND=8) conf_r(larr(3))

       CHARACTER(LEN=larr(12)) dir_raw
!!-----
!! LOCAL VARIABLES
!!-----

       INTEGER(KIND=4) i, j, k, l, m, nid
       INTEGER(KIND=4) n_gal, n_dm, n_aper, n_thread, n_thread2, n_snap, bsize
       REAL(KIND=8) dmp_mass, D2, time(10), pos(3)
       REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp_dbl
       CHARACTER(LEN=100) domnum, fdum, snum, fname
       INTEGER(KIND=4) rd_dlist(larr(2))

       INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: orgind
       !INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: recind
       TYPE(nodetype), DIMENSION(:), ALLOCATABLE :: root
       TYPE(infotype) info
       !TYPE(dat), DIMENSION(:), ALLOCATABLE :: part
       !INTEGER(KIND=4) n_leaf

       TYPE conf_data
         REAL(KIND=8) c_nall, c_n
         REAL(KIND=8) c_mall, c_m
       END TYPE conf_data
       TYPE(conf_data) cdata(larr(3))

       n_gal    = larr(1)
       n_dm    = larr(2)
       n_aper   = larr(3)
       n_thread = larr(4)
       n_thread2= larr(5)
       bsize    = larr(6)
       n_snap   = larr(11)
       
       dmp_mass = darr(1)

       CALL OMP_SET_NUM_THREADS(n_thread)
       !!-----
       !! GET TREE
       !!-----
       IF(ALLOCATED(orgind)) DEALLOCATE(orgind)
       ALLOCATE(orgind(1:n_dm))
       !$OMP PARALLEL DO schedule(static) default(shared)
       DO i=1, n_dm
         orgind(i) = i
       ENDDO
       !$OMP END PARALLEL DO

       time(1) = omp_get_wtime()
       IF(ALLOCATED(root)) DEALLOCATE(root)


       info%bsize = bsize
       info%ndim = 3
       info%npart = n_dm
       info%dtype = 0
       info%vtype = 0
       info%n_thread   = n_thread
       info%dmax_tag   = 1
       !info%omp_tag  = -1

       root = js_kdtree_mktree(xp, mp, orgind, info)
       time(2)  = omp_get_wtime()
       !!-----
       !! MAIN LOOP
       !!-----
       !$OMP PARALLEL DO default(shared) &
       !$OMP & private(cdata, nid, l)
       DO i=1, n_gal
         !! Initialize
         DO l=1, n_aper
           cdata(l)%c_nall = 0.
           cdata(l)%c_mall = 0.
           cdata(l)%c_n = 0.
           cdata(l)%c_m = 0.
         ENDDO
         nid = 1

         !! Walk Tree
         CALL get_contam_walktree(xc(i), yc(i), zc(i), xp, mp, root, nid, aperture(i,:), cdata, n_aper, dmp_mass)
         DO l=1, n_aper
           IF(cdata(l)%c_nall .EQ.0) THEN
             conf_n(i,l) = 0.
           ELSE
             conf_n(i,l) = cdata(l)%c_n / cdata(l)%c_nall
           ENDIF

           IF(cdata(l)%c_mall .EQ.0) THEN
             conf_m(i,l) = 0.
           ELSE
             conf_m(i,l) = cdata(l)%c_m / cdata(l)%c_mall
           ENDIF
         ENDDO

       ENDDO
       !$OMP END PARALLEL DO
       time(3)  = omp_get_wtime()

       PRINT *, time(2)-time(1), ' /  ', time(3)-time(2)
       DEALLOCATE(orgind, root)
       RETURN

CONTAINS
       !! WALK TREE
       RECURSIVE SUBROUTINE get_contam_walktree(x0, y0, z0, xp, mm, root, nid, conf_r, cdata, n_aper, dmp_mass)
       IMPLICIT NONE
       REAL(KIND=8) x0, y0, z0, box(8,3), gpos(3), ppos(3)
       REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: xp
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: mm
       !TYPE(dat), DIMENSION(:), INTENT(IN) :: part
       REAL(KIND=8) dmp_mass
       TYPE(nodetype), DIMENSION(:), INTENT(IN) :: root
       TYPE(nodetype) node
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: conf_r
       TYPE(conf_data) cdata(n_aper)

       INTEGER(KIND=4) n_aper, i, j, k, l, m, nid
       REAL(KIND=8) max_rr, d2, gd2

       !PRINT *, nid, root(nid)%bstart, root(nid)%bend, root(nid)%leaf, root(nid)%dmax
       !! Initialize
       max_rr = conf_r(n_aper)
       gpos(1) = x0; gpos(2) = y0; gpos(3) = z0
       node = root(nid)
       IF(js_kdtree_nodeskip(node, gpos, max_rr, 3)) RETURN
       IF(node%leaf .GT. 0) THEN

         !IF(js_kdtree_nodeskip(node, gpos, max_rr, 3)) RETURN

         DO i=node%bstart, node%bend
           ppos(1) = xp(i,1); ppos(2) = xp(i,2); ppos(3) = xp(i,3)
           gd2 = js_d3d(gpos, ppos)
           DO j=1, n_aper
             IF(gd2 .GT. conf_r(j)) CYCLE
             cdata(j)%c_nall = cdata(j)%c_nall + 1.
             cdata(j)%c_mall = cdata(j)%c_mall + mm(i)

             
             IF(mm(i) .GT. 1.1*dmp_mass) THEN
               cdata(j)%c_n = cdata(j)%c_n + 1.
               cdata(j)%c_m = cdata(j)%c_m + mm(i)
             ENDIF
           ENDDO
         ENDDO
         RETURN
       ELSE

         nid = node%left
         CALL get_contam_walktree(x0, y0, z0, xp, mm, root, nid, conf_r, cdata, n_aper, dmp_mass)

         nid = node%right
         CALL get_contam_walktree(x0, y0, z0, xp, mm, root, nid, conf_r, cdata, n_aper, dmp_mass)
       ENDIF


       RETURN

       END SUBROUTINE get_contam_walktree

       FUNCTION js_d3d(x, y) RESULT(d)
       IMPLICIT NONE
       REAL(KIND=8) x(3), y(3), d
       d = (x(1) - y(1))**2 + (x(2) - y(2))**2 + (x(3) - y(3))**2
       d = SQRT(d)
       RETURN
       END FUNCTION js_d3d

       FUNCTION get_contam_inbox(bnd, gpos, max_r) RESULT(ok)
       IMPLICIT NONE
       REAL(KIND=8) bnd(6,2), gpos(3), max_r
       REAL(KIND=8) box(8,3), d2
       LOGICAL ok
       INTEGER(KIND=4) i
       ok = .False.
       IF(gpos(1) .GE. bnd(1,1) .AND. gpos(1) .LT. bnd(1,2) &
           .AND. gpos(2) .GE. bnd(2,1) .AND. gpos(2) .LT. bnd(2,2) &
           .AND. gpos(3) .GE. bnd(3,1) .AND. gpos(3) .LT. bnd(3,2)) THEN
         ok = .True.
       ELSE
         box(1,1) = bnd(1,1); box(1,2) = bnd(2,1); box(1,3) = bnd(3,1)
         box(2,1) = bnd(1,1); box(2,2) = bnd(2,2); box(2,3) = bnd(3,1)
         box(3,1) = bnd(1,2); box(3,2) = bnd(2,1); box(3,3) = bnd(3,1)
         box(4,1) = bnd(1,2); box(4,2) = bnd(2,2); box(4,3) = bnd(3,1)
  
         box(5,1) = bnd(1,1); box(5,2) = bnd(2,1); box(5,3) = bnd(3,2)
         box(6,1) = bnd(1,1); box(6,2) = bnd(2,2); box(6,3) = bnd(3,2)
         box(7,1) = bnd(1,2); box(7,2) = bnd(2,1); box(7,3) = bnd(3,2)
         box(8,1) = bnd(1,2); box(8,2) = bnd(2,2); box(8,3) = bnd(3,2)
         d2 = js_d3d(box(1,:), gpos)
         DO i=2, 8
           d2 = MIN(d2, js_d3d(box(i,:), gpos))
         ENDDO
         IF(d2 .LT. max_r) ok = .True.
       ENDIF
       RETURN
       END FUNCTION get_contam_inbox

       END SUBROUTINE get_contam

