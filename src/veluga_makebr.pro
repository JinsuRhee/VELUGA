FUNCTION veluga_makebr_init, settings, veluga
	IF settings.makebr_treedir EQ 'des' THEN BEGIN
		treeset	= {$
			N0		: MIN(settings.makebr_snap(0:1)), $
			N1		: MAX(settings.makebr_snap(0:1)), $
			DN		: settings.makebr_snap(2), $
			tag_num	: 'NumDesc', $
			tag_off	: 'DescOffsets', $
			tag_result	: 'Descendants', $
			tag_npart	: 'DescNpart', $
			tag_merit	: 'Merits', $
			tag_nlink	: 'Nsteps_search_new_links' $
			}
	ENDIF ELSE IF settings.makebr_treedir EQ 'prg' THEN BEGIN
		treeset	= {$
			N0		:  MAX(settings.makebr_snap(0:1)), $
			N1		:  MIN(settings.makebr_snap(0:1)), $
			DN		: -settings.makebr_snap(2), $
			tag_num	: 'NumProgen', $
			tag_off	: 'ProgenOffsets', $
			tag_result	: 'Progenitors', $
			tag_npart	: 'ProgenNpart', $
			tag_merit	: 'Merits', $
			tag_nlink	: 'Nsteps_search_new_links' $
			}
	ENDIF ELSE BEGIN
		veluga->ppout, 'Incorrect Tree direction: settings.makebr_treedir'
		STOP
	ENDELSE

	RETURN, treeset
END

;;-----
;; ALLOCATE TREE OBJECT
;;-----
FUNCTION veluga_makebr_allo, settings, veluga, treeset, max_ngal
	snaparr	= LONARR(ABS(treeset.n1-treeset.n0)+1L)-1L
	galarr	= LONARR(5000L)-1L

	dumstr	= {ID:snaparr, snap:snaparr, p_snap:snaparr, p_id:snaparr, p_merit:DOUBLE(snaparr), $
		d_snap:snaparr, d_id:snaparr, endind:-1L, $
		m_id:galarr, m_snap:galarr, m_merit:DOUBLE(galarr), m_BID:galarr, $
		numprog:1L}
	RETURN, REPLICATE(dumstr, max_ngal)
END

;;-----
;; CHANGE THE TREEFROG OUTPUT FILENAME
;;-----
PRO veluga_makebr_tfcname, settings, veluga

	dir     = settings.dir_tree
        flist   = dir + '/tfout/tree.snaplist.txt'
        slist   = LONARR(FILE_LINES(flist))
        OPENR, 10, flist
        FOR i=0L, FILE_LINES(flist)-1L DO BEGIN
                str     = ' '
                READF, 10, str
                ind     = STRPOS(str,'snap_')
                str     = STRMID(str,ind+5,4)
                slist(i)= LONG(str)
        ENDFOR
        CLOSE, 10

	tfile   = FILE_SEARCH(dir + '/tfout/tree.snapshot*.tree')
	IF N_ELEMENTS(tfile) NE N_ELEMENTS(slist) THEN BEGIN
		veluga->errorout, 'Wrong number of tree results: tree.snaplist.txt does not match with the tree output'
		STOP
	ENDIF

	tnumber = LONARR(N_ELEMENTS(tfile))
        FOR i=0L, N_ELEMENTS(tfile)-1L DO BEGIN
                str     = tfile(i)
                str     = STRSPLIT(str, '_' ,/extract)
                str     = str(-1)
                str     = STRSPLIT(str, '.', /extract)
                str     = str(0)
                tnumber(i)      = LONG(str)
        ENDFOR
        cut     = SORT(tnumber)
        tfile   = tfile(cut)

        FOR i=0L, N_ELEMENTS(tfile)-1L DO BEGIN
                orgname = tfile(i)
                newname = 'tree.snapshot_' + STRING(slist(i),format='(I4.4)') + 'VELOCIraptor'
                IF orgname EQ newname THEN CONTINUE
                tmp     = 'mv ' + $
                        orgname + ' ' + dir + '/tfout/' + newname
                SPAWN, tmp
        ENDFOR

        FOR i=0L, N_ELEMENTS(tfile)-1L DO BEGIN
                orgname = 'tree.snapshot_' + STRING(slist(i),format='(I4.4)') + 'VELOCIraptor'
                newname = orgname + '.tree'
                tmp     = 'mv ' + dir + '/tfout/' + orgname + ' ' + dir + '/tfout/' + newname
                SPAWN, tmp
        ENDFOR

END

;;-----
;; FIND adjacent snapshot
;;-----
FUNCTION veluga_makebr_findnextsnap, settings, snap_curr

	snap_next       = snap_curr
        goout   = -1L
        REPEAT BEGIN
                IF settings.makebr_treedir EQ 'des' THEN snap_next ++
                IF settings.makebr_treedir EQ 'prg' THEN snap_next --
                dumfname        = settings.dir_tree + '/tfout/tree.snapshot_' + $
                        STRING(snap_next,format='(I4.4)') + 'VELOCIraptor.tree'
                IF STRLEN(FILE_SEARCH(dumfname)) GE 5L THEN $
                        goout   = 1L
		
                IF snap_next LT settings.makebr_snap(0) OR $
                        snap_next GT settings.makebr_snap(1) THEN BEGIN
                        goout           = 1L
                        snap_next       = -1L
                ENDIF
        ENDREP UNTIL goout GE 1L
        RETURN, snap_next
END

;;-----
;; Link
;;-----
PRO veluga_makebr_link, settings, tree, gind, evoldum, snap_curr, snap_next, $
        dum_id, dum_mer, t_id

        tree(gind).endind ++
        evoldum.id(gind)        = dum_id
        evoldum.snap(gind)      = snap_next
        evoldum.merit(gind)     = dum_mer

        tree(gind).ID(tree(gind).endind)        = t_id
        tree(gind).snap(tree(gind).endind)      = snap_curr
        IF snap_curr EQ tree(gind).snap(tree(gind).endind-1L) THEN STOP ;;123123
        IF settings.makebr_treedir EQ 'des' THEN BEGIN
                tree(gind).p_snap(tree(gind).endind+1L) = snap_curr
                tree(gind).p_ID(tree(gind).endind+1L)           = t_ID
                tree(gind).p_merit(tree(gind).endind+1L)        = dum_mer
                tree(gind).d_snap(tree(gind).endind)            = snap_next
                tree(gind).d_ID(tree(gind).endind)              = dum_id
        ENDIF ELSE IF settings.makebr_treedir EQ 'prg' THEN BEGIN
                tree(gind).p_snap(tree(gind).endind+1)          = snap_next
                tree(gind).p_ID(tree(gind).endind+1)            = dum_id
                tree(gind).p_meirt(tree(gind).endind+1)         = dum_mer
                tree(gind).d_snap(tree(gind).endind)            = snap_curr
                tree(gind).d_ID(tree(gind).endind)              = t_ID
        ENDIF
END

;;-----
;; Link Progenitor
;;-----
PRO veluga_makebr_proglink, settings, tree, ind0, endind, n_comp, p_id, snap_curr, merit
	tree(ind0).numprog ++
        nn      = tree(ind0).numprog-2L
        tree(ind0).m_id(nn)     = p_id
        tree(ind0).m_snap(nn)   = snap_curr
        tree(ind0).m_merit(nn)  = merit
        tree(ind0).m_BID(nn)    = n_comp
END

;;-----
;; Finish Branch
;;-----
PRO veluga_makebr_finishbranch, settings, tree, complete_tree, n_comp, ind, stat
        a       = tree(ind)
        nn      = a.endind
        nn2     = (a.numprog-2L) > 0L
        b       = {ID:a.ID(0L:nn), snap:a.snap(0L:nn), stat:stat, $
                p_snap:a.p_snap(0L:nn), p_id:a.p_id(0L:nn), p_merit:a.p_merit(0L:nn), $
                m_ID:a.m_ID(0:nn2), m_merit:a.m_merit(0:nn2), m_bid:a.m_bid(0L:nn2), m_snap:a.m_snap(0L:nn2), $
                d_snap:a.d_snap(0L:nn), d_id:a.d_id(0L:nn,*), endind:a.endind, numprog:a.numprog}
        complete_tree(n_comp)   = PTR_NEW(b,/no_copy)
        n_comp  ++
END

;;-----
;; CLEAR TREE
;;-----
PRO veluga_makebr_clearbranch, tree, evoldum, ind, gind

        evoldum.ID(ind)         = evoldum.ID(gind)
        evoldum.snap(ind)       = evoldum.snap(gind)
        evoldum.merit(ind)      = evoldum.merit(gind)
        evoldum.ID(gind)        = -1L
        evoldum.snap(gind)      = -1L
        evoldum.merit(gind)     = -1.0d

        tree(ind)       = tree(gind)
        tree(gind).ID           = -1L
        tree(gind).snap         = -1L
        tree(gind).p_snap       = -1L
        tree(gind).p_id         = -1L
        tree(gind).p_merit      = -1.d
        tree(gind).d_snap       = -1L
        tree(gind).d_id         = -1L
        tree(gind).endind       = -1L
        tree(gind).m_id         = -1L
        tree(gind).m_snap       = -1L
        tree(gind).m_merit      = -1.d
        tree(gind).m_bid        = -1L
        tree(gind).numprog      = 1L

        gind --
END

;;----- REALLOCATE
PRO veluga_makebr_reallocate_t, tree, gind, evoldum, maxgind
        maxind  = N_ELEMENTS(tree) + maxgind
        ;; tree
        tree2   = REPLICATE(tree(0), maxind)
        tree2(0L:gind)  = tree(0L:gind)
        tree    = tree2

        ;; evoldum
        evoldum2        = {ID:LONARR(maxind)-1L, snap:LONARR(maxind)-1L, merit:DBLARR(maxind)-1.d}
        evoldum2.ID(0L:gind)    = evoldum.ID(0L:gind)
        evoldum2.snap(0L:gind)  = evoldum.snap(0L:gind)
        evoldum2.merit(0L:gind) = evoldum.merit(0L:gind)
        evoldum = evoldum2
END
PRO veluga_makebr_reallocate_ct, complete_tree, n_comp, maxgind
        maxind  = N_ELEMENTS(complete_tree) + maxgind
        ;; complete_tree
        complete_tree2  = PTRARR(maxind)
        complete_tree2(0L:n_comp-1L)    = complete_tree(0L:n_comp-1L)
        complete_tree   = complete_tree2
END

;;-----
;; READ HDF5
;;-----
FUNCTION veluga_makebr_rdhdf5, fname, treeset

        fid     = H5F_OPEN(fname)

        did     = H5D_OPEN(fid, treeset.tag_num)
        num = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_off)
        off = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_result)
        res = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_merit)
        mer = H5D_READ(did) & H5D_CLOSE, did

        did     = H5D_OPEN(fid, treeset.tag_npart)
        npart = H5D_READ(did) & H5D_CLOSE, did

        did     = H5A_OPEN_NAME(fid, treeset.tag_nlink)
        nlink = H5A_READ(did) & H5A_CLOSE, did

        did     = H5D_OPEN(fid, 'ID')
        id = H5D_READ(did) & H5D_CLOSE, did

        H5F_CLOSE, fid

        dum     = {num:num, off:off, res:res, merit:mer, npart:npart, nlink:nlink, id:id}
        RETURN, dum
END

;;-----
;; Arrange the finished branches
;;-----
PRO veluga_makebr_lastsnap, settings, veluga, tree, complete_tree, n_comp, evoldum, gind, g_curr, snap_curr

	IF settings.makebr_treedir EQ 'prog' THEN BEGIN
		veluga->errorout, 'progenitor direction is not implemented yet'
                STOP
        ENDIF

	;;-----
	;; Scan by IDs of the last snapshot
	;;-----
	FOR i=0L, N_ELEMENTS(g_curr)-1L DO BEGIN
		cut_exist	= WHERE(evoldum.snap EQ snap_curr AND $
			evoldum.id EQ g_curr(i).ID, ncut_exist)

		IF ncut_exist EQ 1L THEN BEGIN  ;; Finish Single Branch
                        tree(cut_exist).endind ++
                        tree(cut_exist).ID(tree(cut_exist).endind)      = g_curr(i).ID
                        tree(cut_exist).snap(tree(cut_exist).endind)    = snap_curr
		ENDIF ELSE IF ncut_exist GE 2L THEN BEGIN
			;; Merger happened at the last snapshot
			;;	Select branch with higher merit as the main branch
			maxmerit        = MAX(evoldum.merit(cut_exist))
                        ind0            = cut_exist(WHERE(evoldum.merit(cut_exist) EQ maxmerit))
                        endind          = tree(ind0).endind + 1
			FOR li=0L, N_ELEMENTS(cut_exist)-1L DO BEGIN
                                ind     = cut_exist(li)
                                IF evoldum.merit(ind) EQ maxmerit THEN BEGIN
                                        ;; Primary Link
                                        tree(ind).endind ++
                                        tree(ind).ID(tree(ind).endind)  = g_curr(i).ID
                                        tree(ind).snap(tree(ind).endind)= snap_curr
                                ENDIF ELSE IF ABS(tree(ind).snap(0) - snap_curr) GE 10L AND $
                                        tree(ind).snap(0) NE -1L THEN BEGIN
                                        ;; Secondary Link
                                        ;;      Merit is low but has a reasonable tree length
                                        veluga_makebr_proglink, settings, tree, ind0, endind, n_comp, $
                                                tree(ind).id(tree(ind).endind), $
                                                tree(ind).snap(tree(ind).endind), $
                                                evoldum.merit(ind)
                                        veluga_makebr_finishbranch, settings, tree, complete_tree, n_comp, ind, 'sub'
                                        veluga_makebr_clearbranch, tree, evoldum, ind, gind
                                ENDIF ELSE BEGIN
                                        veluga_makebr_clearbranch, tree, evoldum, ind, gind
                                ENDELSE

                        ENDFOR
		ENDIF
	ENDFOR

	tree	= tree(0L:gind)

END

;;-----
;; MATCHING
;;-----
PRO veluga_makebr_match, settings, veluga, treelog, tree, complete_tree, n_comp, t_curr, g_curr, g_next, gind, evoldum, snap_curr

	mlimit  = settings.makebr_meritlimit

	;; loop by galaxy
	FOR i=0L, N_ELEMENTS(t_curr.id)-1L DO BEGIN
		IF t_curr.num(i) EQ 0L THEN CONTINUE

		treelog.n_all ++

		ind1    = t_curr.off(i)
                ind2    = t_curr.off(i) + t_curr.num(i) - 1L
                dum_id  = t_curr.res(ind1:ind2)
                dum_part= t_curr.npart(ind1:ind2)
                dum_mer = t_curr.merit(ind1:ind2)

		;;----- Matching by one with most largest Merit
                cut     = MIN(WHERE(dum_mer EQ MAX(dum_mer)))
                cut2    = WHERE(g_next.ID EQ dum_id(cut) AND g_next.npart EQ dum_part(cut), ncut)

		IF ncut EQ 0L THEN CONTINUE	;; no link

		IF ncut GE 2L THEN BEGIN        ;;----- For a galaxy that is not changed in its particle number over the multiple snapshots
                        nochagal_snap   = g_next(cut2).snapnum
                        cut_nochagal    = WHERE(nochagal_snap EQ MIN(nochagal_snap))
                        cut2    = cut2(cut_nochagal)
                ENDIF

		snap_next       = g_next(cut2).snapnum


		;;----- Link

		cut_exist       = WHERE(evoldum.snap EQ snap_curr AND $
                        evoldum.id EQ t_curr.ID(i), ncut_exist) ;; branch exist?
		IF ncut_exist EQ 0L THEN BEGIN  ;; NEW MERGER TREE
                        gind ++
                        veluga_makebr_link, settings, tree, gind, evoldum, snap_curr, snap_next, $
                                dum_id(cut), dum_mer(cut), t_curr.ID(i)

                        treelog.n_new ++
                ENDIF ELSE BEGIN		;; BRANCH EXIST

			;; IF the second merit is too low, just merged
                        maxmerit        = MAX(evoldum.merit(cut_exist))
                        ind0            = MIN(cut_exist(WHERE(evoldum.merit(cut_exist) EQ maxmerit)))
                        endind          = tree(ind0).endind+1

			FOR li=0L, N_ELEMENTS(cut_exist)-1L DO BEGIN
				ind     = cut_exist(li)
                                IF ind EQ ind0 THEN BEGIN;evoldum.merit(ind) EQ maxmerit THEN BEGIN
                                        ;; Primary Link
                                        veluga_makebr_link, settings, tree, ind, $
                                                evoldum, snap_curr, snap_next, $
                                                dum_id(cut), dum_mer(cut), t_curr.ID(i)
                                        treelog.n_link ++
				ENDIF ELSE IF ABS(tree(ind).snap(0)-snap_curr) GE 10L AND $
                                        tree(ind).snap(0) NE -1L THEN BEGIN
                                        ;; Secondary Link
                                        ;;      Merit is low but has a reasonable tree length
                                        veluga_makebr_proglink, settings, tree, ind0, endind, n_comp, $
                                                tree(ind).id(tree(ind).endind), $
                                                tree(ind).snap(tree(ind).endind), $
                                                evoldum.merit(ind)
                                        veluga_makebr_finishbranch, settings, tree, complete_tree, n_comp, ind, 'sub'

                                        store_ind       = WHERE(cut_exist EQ gind, nstore)
                                        veluga_makebr_clearbranch, tree, evoldum, ind, gind
                                        treelog.n_link2 ++
                                        IF nstore GE 1L THEN cut_exist(store_ind) = ind
				ENDIF ELSE BEGIN
                                        ;; Too short branch or too low merit
                                        ;;      Clear it
                                        store_ind       = WHERE(cut_exist EQ gind, nstore)

                                        veluga_makebr_clearbranch, tree, evoldum, ind, gind
                                        treelog.n_broken ++

                                        IF nstore GE 1L THEN cut_exist(store_ind) = ind
                                ENDELSE
                        ENDFOR
		ENDELSE
	ENDFOR

END
;;-----
;; REMOVE FROM THE MAIN BRANCH
;;-----
PRO veluga_makebr_remove, settings, evoldum, tree, gind, complete_tree, n_comp, snap_curr, nlink

        gind2   = gind
        FOR i=0L, gind2 DO BEGIN
                IF i GT gind THEN BREAK
                IF ABS(tree(i).snap(tree(i).endind)-snap_curr) GE nlink+1L THEN BEGIN
                        IF ABS(tree(i).snap(0)-snap_curr) GE 10L THEN BEGIN
                                veluga_makebr_finishbranch, settings, tree, complete_tree, n_comp, i, 'main'
                                veluga_makebr_clearbranch, tree, evoldum, i, gind
                        ENDIF ELSE BEGIN
                                veluga_makebr_clearbranch, tree, evoldum, i, gind
                        ENDELSE
                ENDIF
        ENDFOR
END

;;----- GEN KEY
FUNCTION veluga_makebr_genkey, settings, tree

        MAX_snap        = 200L
        MAX_ID          = 10000L

        genkey_redo:
        tree_key        = LONARR(MAX_snap + settings.makebr_bidkey*MAX_ID) - 1L

        n_tree          = N_ELEMENTS(tree)
        FOR i=0L, n_tree-1L DO BEGIN
                tmp             = *tree(i)
                s               = tmp.snap
                id              = tmp.id
                ind     = s + settings.makebr_bidkey*id

                IF MAX(id) GT MAX_ID THEN BEGIN
                        MAX_ID          = MAX(id) + 1L
                        GOTO, genkey_redo
                ENDIF
                IF MAX(s) GT MAX_snap THEN BEGIN
                        MAX_snap        = MAX(s) + 1L
                        GOTO, genkey_redo
                ENDIF

                IF MAX(s) GT settings.makebr_bidkey THEN BEGIN
                        power   = LOGN(ALOG10(MAX(s))) + 1.d
                        settings.P_makebr_bidkey = LONG(10.d^power)
                        GOTO, genkey_redo
                ENDIF


                tree_key(ind)   = i
        ENDFOR
        RETURN, tree_key
END

;;-----
;; MAIN
;;-----
PRO veluga_makebr, header, num_thread=num_thread, horg=horg

	;;-----
	;; CALL OBJECT
	;;-----
	veluga	= OBJ_NEW('veluga', header, num_thread=num_thread)

	settings 	= veluga->getheader()
	settings 	= CREATE_STRUCT(settings, 'horg', horg)

	IF horg EQ 'g' THEN BEGIN
		dir_tree	= settings.dir_catalog + 'Galaxy/tree'
	ENDIF ELSE BEGIN
		dir_tree	= settings.dir_catalog + 'Halo/tree'
	ENDELSE
	settings	= CREATE_STRUCT(settings, 'dir_tree', dir_tree)
	
	treeset 	= veluga_makebr_init(settings, veluga)

	;;veluga_makebr_initcompile, settings, veluga

	;;-----
	;; ALLOCATE
	;;-----
	max_ngal	= 10000L
		;; size of branch array (automatically reallocate when N>max_ngal) updated?
	tree	= veluga_makebr_allo(settings, veluga, treeset, max_ngal)
		;; array for branch

	complete_tree	= PTRARR(max_ngal)
		;; complete branch array

	evoldum	= {ID:LONARR(max_ngal)-1L, snap:LONARR(max_ngal)-1L, merit:DBLARR(max_ngal)-1.d}

	;;-----
	;; CHANGE THE TREEFROG OUTPUT FILENAME
	;;-----
	veluga_makebr_tfcname, settings

	;;-----
	;; MAIN LOOP
	;;-----
	gind	= -1L
	n_comp	= 0L
	treelog	= {n_new:0L, n_link:0L, n_link2:0L, n_link3:0L, n_broken:0L, n_all:0L}

	FOR i=treeset.n0, treeset.n1, treeset.dn DO BEGIN
		;; SAVE PART?
		;;----- SNAPSHOT CHECK (skip if this snapshot is empty)
		dumfname        = settings.dir_catalog
                IF settings.horg EQ 'h' THEN dumfname += 'Halo/'
                IF settings.horg EQ 'g' THEN dumfname += 'Galaxy/'
                dumfname += 'snap_' + STRING(i,format='(I4.4)')

		IF STRLEN(FILE_SEARCH(dumfname)) LE 5L THEN CONTINUE

		snap_curr	= i
		snap_next	= veluga_makebr_findnextsnap(settings, snap_curr)

		;;----- LAST SNAPSHOT
		IF i EQ treeset.n1 OR snap_next EQ -1L THEN BEGIN
			g_curr	= veluga->r_gal(snap_curr, -1L, Gprop=['ID', 'npart'], horg=settings.horg)

			veluga_makebr_lastsnap, settings, veluga, tree, compelete_tree, n_comp, evoldum, gind, g_curr, snap_curr


			FOR j=0L, gind DO BEGIN
				IF N_ELEMENTS(complete_tree) - n_comp LT 2000L THEN $
                                        veluga_makebr_reallocate_ct, complete_tree, n_comp, maxgind

                                IF ABS(tree(j).snap(0)-snap_curr) GE 10L THEN $
                                        veluga_makebr_finishbranch, settings, tree, complete_tree, n_comp, j, 'main'
			ENDFOR

			complete_tree   = complete_tree(0L:n_comp-1L)
			BREAK
		ENDIF

		;;----- LOAD REQUIRED DATA
		treefname	= settings.dir_tree + '/tfout/tree.snapshot_' + STRING(snap_curr, format='(I4.4)') + 'VELOCIraptor.tree'

		t_curr	= veluga_makebr_rdhdf5(treefname, treeset)
		g_curr	= veluga->r_gal(snap_curr, -1L, Gprop=['ID', 'npart'], horg=settings.horg)
		g_next	= veluga->r_gal(snap_next, -1L, Gprop=['ID', 'npart'], horg=settings.horg)


		IF t_curr.nlink GE 2L THEN BEGIN
			FOR i2=0L, t_curr.nlink-2L DO BEGIN
				snap_next = veluga_makebr_findnextsnap(settings, snap_next)
				IF snap_next EQ -1L THEN CONTINUE
				g_dum	= veluga->r_gal(snap_next, -1L, Gprop=['ID', 'npart'], horg=settings.horg)
				g_next	= [g_next, g_dum]
			ENDFOR
		ENDIF


		;;----- MATCHING
		veluga_makebr_match, settings, veluga, treelog, tree, complete_tree, n_comp, t_curr, g_curr, g_next, gind, evoldum, snap_curr

		;;----- REMOVE BRANCH
		veluga_makebr_remove, settings, evoldum, tree, gind, complete_tree, n_comp, snap_curr, t_curr.nlink

		;;----- TREELOG
		PRINT, i, ' / ', treeset.N1, ' Using n_step = ', t_curr.nlink
                PRINT, '        ALL : ' +  STRTRIM(treelog.n_all,2) + $
                        ' / NEW :' + STRTRIM(treelog.n_new,2) +  $
                        ' / LINK : ' + STRTRIM(treelog.n_link,2) + $
                        ;' ' + STRTRIM(treelog.n_link2,2) + $
                        ;' ' + STRTRIM(treelog.n_link3,2) + $
                        ' / Broken : ' + STRTRIM(treelog.n_broken,2) + $
                        ' / NGal : ' + STRTRIM(gind+n_comp,2) + $
                        ' / Gind : ' + STRTRIM(gind,2)

                FOR ii=0L, N_TAGS(treelog)-1L DO treelog.(ii) = 0L
	ENDFOR

	;; GENERATE KEY
	tree_key        = veluga_makebr_genkey(settings, complete_tree)
	tree_key(0)	= settings.makebr_bidkey

	;; KEY ASSIGNMENT CHECK
	FOR i=0L, N_ELEMENTS(complete_tree)-1L DO BEGIN
                tt      = *complete_tree(i)
                IF TYPENAME(tt) EQ 'UNDEFINED' THEN CONTINUE
                keyval  = tt.snap + tree_key(0)*tt.id
                keyind  = tree_key(keyval)
                void    = WHERE(keyind NE i, nn)
                IF nn GE 1L THEN BEGIN
                        veluga->errorout, 'WRONG ASSIGNMENT OF KEY'
                        STOP
                ENDIF
        ENDFOR

	SAVE, filename=settings.dir_tree + 'tree.sav', complete_tree, tree_key

END
