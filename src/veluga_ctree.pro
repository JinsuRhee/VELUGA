;;-----
;; ALLOCATE
;;------
FUNCTION veluga_ctree_allocate, settings, complete_tree, tree_key, nn

	tmp     = {ID0:0L, snap0:0L, $      ;; ID & snap at the final snapshot
        ID:0L, snap:0L, $       		;; End of id and snap at the current tree
        pos:DBLARR(3), vel:DBLARR(3), $ ;; position & velocity
        stat:'T', $             		;; T (Tree exists), C (To be connected), B (Broken)
        detstat:-1L, $          		;; -1 (End is not specified yet) 1 (specified)
        p_list:PTR_NEW(0.d), $  		;; Particle list at the ending point
        n_ptcl:-1L, $           		;; # of p_list
        list:REPLICATE({merit:0.d, id:-1L, snap:-1L}, settings.ctree_n_search), $	;; Meirt list
        list_n:0L, $
        etc:'etc'}

    tmp		= REPLICATE(tmp, nn)

    RETURN, tmp
END

;;-----
;; GET TREE
;;-----
FUNCTION veluga_ctree_gtree, ctree, key, snap0, id0

	keyval 	= snap0 + id0*key(0)
	ind 	= key(keyval)
	IF ind EQ -1L THEN RETURN, 1L
	RETURN, *ctree(ind)
END
;;-----
;; REALLOCATE
;;-----
FUNCTION veluga_ctree_reallocate, data, nn

	data2	= REPLICATE(data(0), nn)
	data2(0L:N_ELEMENTS(data)-1L) 	= data
	data 	= 0.d
	RETURN, data2
END
;;-----
;; Input galaxy data into data
;;-----
PRO veluga_ctree_intputgal, settings, complete_tree, tree_key, data, gal, i0, i1

	veluga 	= settings.veluga

	IF i1-i0+1L NE N_ELEMENTS(gal) THEN BEGIN
		veluga->errorout, 'WRONG NUMBER OF GALAXIES WITH INDICES'
		STOP
	ENDIF

	IF i1 GT N_ELEMENTS(data)-1L THEN data = veluga_ctree_reallocate(data, i1 + N_ELEMENTS(gal))
	
	data(i0:i1).id0 	= gal.id
	data(i0:i1).snap0 	= gal.snapnum
	FOR i=0L, N_ELEMENTS(gal)-1L DO BEGIN
		tree 	= veluga_ctree_gtree(complete_tree, tree_key, gal(i).snapnum, gal(i).id)

		IF TYPENAME(tree) EQ 'LONG' THEN BEGIN
			data(i0+i).snap 	= gal(i).snapnum
			data(i0+i).id 		= gal(i).id
		ENDIF ELSE BEGIN
			data(i0+i).snap 	= tree.snap(0)
			data(i0+i).id 		= tree.id(0)
		ENDELSE
	ENDFOR
END

;;-----
;; Get snapshot info
;;-----
FUNCTION veluga_ctree_getsnapinfo, settings, veluga
	;; Get snapshot list

	file 	= settings.dir_tree + '/tfout/tree.snapshot_*VELOCIraptor.tree'
	flist	= FILE_SEARCH(file)
	slist	= LONARR(N_ELEMENTS(flist))

	FOR i=0L, N_ELEMENTS(flist) -1L DO BEGIN
		tmp     = flist(i)
        tmp     = (STRSPLIT(tmp, '/', /extract))[-1]
        i0      = STRPOS(tmp, '_')
        i1      = STRPOS(tmp, 'VELO')
        tmp0    = STRMID(tmp, i0+1L, i1-i0-1L)
        slist(i)= LONG(tmp0)
	ENDFOR

	sinfo 	= REPLICATE({snap:-1L, aexp:0.d, unit_l:0.d, unit_t:0.d, age:0.d, kpc:0.d}, MAX(slist)+1L)
	;; to be index=snapshot
	sinfo(slist).snap 	= slist

	FOR i=0L, N_ELEMENTS(sinfo)-1L DO BEGIN
		IF sinfo(i).snap LT 0L THEN CONTINUE
		
		info 	= veluga->g_info(i)
		sinfo(i).aexp 	= info.aexp
		sinfo(i).unit_l	= info.unit_l
		sinfo(i).unit_t	= info.unit_t
		sinfo(i).age 	= veluga->g_ztoage(1./sinfo(i).aexp-1.d, info=info)
		sinfo(i).kpc 	= info.cgs.kpc
	ENDFOR

	RETURN, {slist:slist, sinfo:sinfo}
END

;;-----
;; TREE CLASSIFICATION
;;-----
PRO veluga_ctree_classify, settings, data, snap0, number

	slist 	= settings.slist

	number 	= {T:0L, C:0L, B:0L}
	ind0    = (WHERE(slist EQ snap0))[0]
	FOR i=0L, N_ELEMENTS(data)-1L DO BEGIN
		;; Tree exists

		IF data(i).snap LE snap0 THEN BEGIN
			data(i).stat 	= 'T'
			number.T ++
			CONTINUE
		ENDIF

		ind     = (WHERE(slist EQ data(i).snap))[0]

		;; For altered Tree
        IF data(i).list(0).snap GT 0L THEN BEGIN
                ind1    = (WHERE(slist EQ data(i).list(0).snap))[0]
        ENDIF ELSE BEGIN
                ind1    = ind0 + settings.ctree_n_search*2L
        ENDELSE

        IF ind - ind0 LE settings.ctree_n_search OR ind1 - ind0 LT settings.ctree_n_search THEN BEGIN
                data(i).stat    = 'C'
                number.C        ++
        ENDIF ELSE BEGIN
                data(i).stat    = 'B'
                number.B        ++
        ENDELSE

	ENDFOR

	PRINT, ''
    PRINT, '                        TREE CLASSIFICATION'
    PRINT, '                                With tree       = ', number.T
    PRINT, '                                To be connected = ', number.C
    PRINT, '                                Broken          = ', number.B
    PRINT, ''
END

;;-----
;; Determine ending point
;;-----
PRO velgua_ctree_detend, settings, data, complete_tree, tree_key
	;snapinfo 	= settings.snapinfo
	cut     = WHERE(data.stat EQ 'C', ncut)

    IF ncut GE 1L THEN BEGIN
        ;FOR i=0L, ncut-1L DO BEGIN
        ;    ind     = cut(i)
        ;    tkey    = data(ind).snap + data(ind).id*tree_key(0)
        ;    tind    = tree_key(tkey)
        ;ENDFOR
        data(cut).detstat       = 1L
    ENDIF
END

;;-----
;; Give the weight to particles
;;		following the formula in Elahi+19b when particles are stored in their binding energy order
;;------
FUNCTION veluga_ctree_getweight, pid

	nn 	= N_ELEMENTS(pid)
	weight	= DINDGEN(nn)+1.d
	weihgt	= REVERSE(weight) / nn
	weight 	/= (0.5772156649d + ALOG(nn*1.d))

	RETURN, weight
END

;;-----
;; Collect particle ID along the branch
;;-----
FUNCTION veluga_ctree_collectpidalongbranch, settings, slist, idlist
	veluga 	= settings.veluga

	;pid 	= veluga->r_pid(slist(0), idlist(0), horg=settings.horg)
	;pweight0= veluga_ctree_getweight(pid)

	;;----- Allocate Large array

	nn 		= (veluga->r_gal(slist(0), idlist(0), horg=settings.horg, GProp=['npart'])).npart

	pid 	= LON64ARR( (10LL*nn*settings.ctree_n_step_N) < 500000000LL )
	pweight = DBLARR( (10LL*nn*settings.ctree_n_step_N) < 500000000LL )

	ind0	= 0L
	loop_n 	= 0L
	loop_ind= 0L

	REPEAT BEGIN
		;;----- IF this branch has a short tree, stop here
		IF loop_ind GE N_ELEMENTS(slist) THEN BREAK

		pid0	= veluga->r_pid(slist(loop_ind), idlist(loop_ind), horg=settings.horg)
		pweight0= veluga_ctree_getweight(pid0)


		ind1 	= ind0 + N_ELEMENTS(pid0)-1L

		;; REALLOCATE IF the large array is full
		IF ind1 GE N_ELEMENTS(pid) THEN BEGIN
			pid 	= [pid, LON64ARR((ind1-ind0)*10L)]
			pweight = [pweight, DBLARR((ind1-ind0)*10L)]			
		ENDIF

		pid(ind0:ind1) 		= pid0
		pweight(ind0:ind1)	= pweight0
		ind0 	= ind1 + 1L

		loop_n ++
		loop_ind 	+= settings.ctree_n_step_dN
	ENDREP UNTIL loop_n GE settings.ctree_n_step_N OR loop_ind GE N_ELEMENTS(slist)

	pid 	= pid(0L:ind0-1LL)
	pweight = pweight(0L:ind0-1LL)

	;;----- REMOVE PARTICLES SHOWING MULTIPLE TIMES
	sortind = SORT(pid)
    pid     = pid(sortind)
    pweight = pweight(sortind)

    uind    = UNIQ(pid)
    pid     = pid(uind)
    pweight = pweight(uind)

    ;;----- LEAVE PARTICLES SHOWN MULTIPLE TIMES, which means they are bound to this galaxy for a long time
    ;;	only 25% of particles which have a long history in this galaxy is used
    ;; MAGIC trick is used here ^^

    numid   = uind - [-1L, uind(0L:-2L)]
	n_step_bw0      = settings.ctree_n_step_N
    REPEAT BEGIN
    	ucut    = WHERE(numid GE n_step_bw0, uncut);tree_set.n_step_bw/denu, uncut)
        n_step_bw0 --
    ENDREP UNTIL DOUBLE(uncut) / DOUBLE(N_ELEMENTS(pid)) GT 0.25d

    pid 	= pid(ucut)
    pweight = pweight(ucut)
    RETURN, {pid:pid, weight:pweight, n_con:n_step_bw0+1L}
END

;;-----
;; Collect particle ID
;;-----
FUNCTION veluga_ctree_collectpid, settings, data, complete_tree, tree_key

	veluga 	= settings.veluga
	sinfo 	= settings.sinfo
	slist 	= settings.slist
	
	;;----- Find galaxies that will be connected and their particles are not collected
	cut     = WHERE(data.stat EQ 'C' AND data.n_ptcl LT 0L, ncut)

	IF ncut GE 1L THEN BEGIN
        FOR i=0L, ncut-1 DO BEGIN
            ind     = cut(i)

            kval    = data(ind).snap + tree_key(0)*data(ind).id
            tind    = tree_key(kval)

            IF tind LT 0L THEN BEGIN
            	;;----- no tree for this galaxy
                pid0    = veluga->r_pid(data(ind).snap, data(ind).id, horg=settings.horg)
				pweight0= veluga_ctree_getweight(pid0) ;; weight of ptcls for merit calculation
			ENDIF ELSE BEGIN
				;;----- This galaxy has a branch. Collect particles along its branch
                tree    = *complete_tree(tind)
                t_cut   = WHERE(tree.snap GE data(ind).snap, t_nn)

                IF t_nn EQ 0L THEN STOP ;; Weird branch. Stop

                ;;----- Snapshot / Galaxy ID list along a branch before its branch end
                t_slist = tree.snap(t_cut)
                t_idlist= tree.id(t_cut)

                cpid    = veluga_ctree_collectpidalongbranch(settings, t_slist, t_idlist)

				pid0    = cpid.pid
                pweight0= cpid.weight
            ENDELSE

			;;----- STORE PARTICLE TO THIS BRANCH ARRAY
			data(ind).n_ptcl        = N_ELEMENTS(pid0)
            data(ind).p_list        = PTR_NEW(pid0, /no_copy)

            ;;----- STORE GALAXY 6D Coordinates (to read the raw particles within a volume)
			gal0    = veluga->r_gal(data(ind).snap, data(ind).id, horg=settings.horg, Gprop=['Xc', 'Yc', 'Zc', 'VXc', 'VYc', 'VZc'])
            data(ind).pos   = [gal0.xc, gal0.yc, gal0.zc]
            data(ind).vel   = [gal0.vxc, gal0.vyc, gal0.vzc]
        ENDFOR
    ENDIF

    ;;
    cut 	= WHERE(data.n_ptcl GE 0L, ncut)
    IF ncut EQ 0L THEN BEGIN
    	veluga->errorout, 'No branch has particles?'
    	STOP
    ENDIF

    ;;----- COLLECT All galaxy particles to make a single array
    ;;	store their host galaxy id simultaneously
    nptcl 	= TOTAL(data(cut).n_ptcl)

    dum     = {pid:LONG64(1), gid:0L}
    pid     = REPLICATE(dum, nptcl)

    i0      = 0L
    FOR i=0L, ncut-1L DO BEGIN
        ind     = cut(i)

        i1      = i0+data(ind).n_ptcl-1L
        IF i1 GE N_ELEMENTS(pid) THEN BEGIN
            pid     = [pid, REPLICATE(pid(0),i1-i0+1L)]
        ENDIF
        pid(i0:i1).pid  = *data(ind).p_list
        pid(i0:i1).gid  = data(ind).id0

        i0      = i1+1L
    ENDFOR
    RETURN, pid(0L:i1)
END

;;-----
;; READ Snapshot particles
;;-----
FUNCTION veluga_ctree_readsnap, settings, data, snap0

	sinfo 	= settings.sinfo
	slist 	= settings.slist
	veluga 	= settings.veluga

	;;-----
	;; FIRST FIND COMOVING VOLUME
	;;-----
	cut     = WHERE(data.stat EQ 'C', ncut)
    cen     = DBLARR(ncut,3)
    rad     = DBLARR(ncut)
    FOR i=0L, ncut-1L DO BEGIN
        ind     = cut(i)
        pos     = data(ind).pos
        vel     = data(ind).vel

        
        unit_l  = sinfo(data(ind).snap).unit_l
        kpc 	= sinfo(data(ind).snap).kpc

        cen(i,*)= pos*kpc/unit_l

        speed   = NORM(vel) / (kpc/1d5) * (365.d * 86400d * 1e9)        ;; [kpc/Gyr]
        rad(i)  = speed * settings.ctree_rfact * $
        	ABS(sinfo(snap0).age - sinfo(data(ind).snap).age) * kpc / unit_l
    ENDFOR

    gal 	= veluga->r_gal(snap0, -1L, horg=settings.horg, GProp=['ID', 'npart', 'Xc', 'Yc', 'Zc'])

    nptcl 	= TOTAL(gal.npart)
    dum 	= {pid:LONG64(0), gid:-1L}
    pid 	= REPLICATE(dum, nptcl)

    i0 	= 0L
    checkind        = LONARR(N_ELEMENTS(gal))-1L
    FOR i=0L, N_ELEMENTS(gal)-1L DO BEGIN
        ;;----- Distance check

        unit_l  = sinfo(snap0).unit_l
        kpc 	= sinfo(snap0).kpc


        pos     = [gal(i).xc, gal(i).yc, gal(i).zc] * kpc / unit_l
        d3d     = veluga->g_d3d(cen(*,0), cen(*,1), cen(*,2), pos) / rad

        IF MIN(d3d) GT 1. THEN CONTINUE
        checkind(i)     = 1L

        pid0    = veluga->r_pid(snap0, gal(i).id, horg=settings.horg)

        npart0  = N_ELEMENTS(pid0)
        i1      = i0+npart0-1L
        IF i1 GE N_ELEMENTS(pid) THEN pid = [pid, REPLICATE(dum,npart0)]

        pid(i0:i1).pid  = pid0
        pid(i0:i1).gid  = gal(i).id
        i0      = i1 + 1L
    ENDFOR

    cut     = WHERE(pid.gid GT 0L, ncut)
    IF ncut EQ 0L THEN BEGIN
        veluga->errorout, 'NO galaxies within this volume. Increase c_tree_rfact'
        STOP
    ENDIF
    pid     = pid(cut)

    ;;-----
    ;; Make Hash
    ;;-----


    ;; split index
    ;;  There is a memory bug in the OMP implementation
    ;;  Since a serial calculation shows a good performance, bug is remained to be fixed
    num_thread      = 1L


    dN      = N_ELEMENTS(pid)/num_thread
    indarr  = LONARR(num_thread, 2L)
    
    FOR i=0L, num_thread-1L DO BEGIN
        indarr(i,0)     = dN*i
        indarr(i,1)     = dN*(i+1L) - 1L
    ENDFOR
    indarr(-1,1)    = N_ELEMENTS(pid)-1L

    ;; hash allocate
    dN      = MAX(indarr(*,1)-indarr(*,0))+1L
    hash            = LONARR(dN,num_thread) - 1L
    hash_next       = LONARR(dN,num_thread) - 1L

    ;; fortran
    ftr_name        = settings.dir_lib + 'src/fortran/get_hash.so'
        larr = LONARR(20) & darr = DBLARR(20)
        larr(0) = N_ELEMENTS(pid)
        larr(1) = dN

        larr(10)= num_thread

    void    = CALL_EXTERNAL(ftr_name, 'get_hash', $
        larr, darr, pid.pid, indarr, hash, hash_next)


    cut     = WHERE(checkind GE 0L)

    RETURN, {pid:pid, n_ptcl:gal(cut).npart, gid:gal(cut).id, hash:hash, hash_next:hash_next, dn:dn}

END

;;-----
;; MERIT CALCULATION
;;-----
PRO veluga_ctree_commerit, settings, data, pid, pid0, c_snap
	pid_g   = pid.pid
    gid_g   = pid.gid

    pid_s   = pid0.pid.pid
    gid_s   = pid0.pid.gid
    hash            = pid0.hash
    hash_next       = pid0.hash_next


    ;; ALLOCATE (ID = IND in fortran)
    merit   = DBLARR(MAX(gid_g)+1L, MAX(gid_s)+1L)
    npart_g = LONARR(MAX(gid_g)+1L)
    npart_s = LONARR(MAX(gid_s)+1L)

    ;; PTCL NUM INPUT
    cut     = WHERE(data.stat EQ 'C', ncut)
    npart_g(data(cut).ID0-1L)       = data(cut).n_ptcl
    npart_s(pid0.gid-1L)            = pid0.n_ptcl

    match_id        = LONARR(N_ELEMENTS(cut)) - 1L
    match_merit     = DBLARR(N_ELEMENTS(cut)) - 1.d

    ftr_name        = settings.dir_lib + 'src/fortran/get_merit2.so'
        larr = LONARR(20) & darr = DBLARR(20)
        larr(0) = N_ELEMENTS(pid_g)
        larr(1) = N_ELEMENTS(pid_s)
        larr(2) = N_ELEMENTS(npart_g)
        larr(3) = N_ELEMENTS(npart_s)
        larr(4) = pid0.dn
        larr(5) = N_ELEMENTS(hash(0,*))
        larr(6) = N_ELEMENTS(cut)

        larr(10)= 1L;tree_set.num_thread

        ;;-----
        ;; Again, a serial calculation has a good performance, but OMP has a bug at the moment
        

    void    = CALL_EXTERNAL(ftr_name, 'get_merit2', $
    	larr, darr, pid_g, gid_g, pid_s, gid_s, $
        hash, hash_next, $
        npart_g, npart_s, merit, $
        match_id, match_merit)

    FOR i=0L, ncut-1L DO BEGIN
        n0      = data(cut(i)).list_n
        data(cut(i)).list(n0).merit     = match_merit(i)
        data(cut(i)).list(n0).id        = match_id(i)
        data(cut(i)).list(n0).snap      = c_snap
        data(cut(i)).list_n ++
    ENDFOR

END
;;-----
;; LINK
;;-----
PRO veluga_ctree_link, settings, data, number, c_snap, complete_tree, tree_key
	number 	= {n_link:0L, n_broken:0L}
	slist0 	= settings.slist

	ind1    = (WHERE(slist0 EQ c_snap))[0]
    ind0    = (WHERE(slist0 EQ settings.ctree_snap(0)))[0]
    snap_int_cut    = MIN([ind1-ind0-1L, settings.ctree_n_search])
    cut     = WHERE(data.stat EQ 'C' AND data.list_n GE snap_int_cut AND data.list_n GT 0L, ncut)
    IF ncut EQ 0L THEN RETURN

    FOR i=0L, ncut-1L DO BEGIN
        ind     = cut(i)
        mlist   = data(ind).list.merit
        idlist  = data(ind).list.id
        slist   = data(ind).list.snap
        IF MAX(mlist) LT settings.ctree_meritlimit THEN BEGIN
            data(ind).stat = 'B'
            veluga_ctree_free, data, ind, -1L, -1L, c_snap
            CONTINUE
        ENDIF
        
        cut2    = (WHERE(mlist EQ MAX(mlist)))[0]

        id_tolink       = idlist(cut2)
        snap_tolink     = slist(cut2)
        merit_tolink    = MAX(mlist)
        
        ;; Check Whether there is tree existed
        kval    = snap_tolink + tree_key(0)*id_tolink
        tind    = tree_key(kval)

        IF tind EQ -1L THEN BEGIN       ;; notree
            veluga_ctree_expandbr, data, ind, complete_tree, tree_key, id_tolink, snap_tolink, merit_tolink

            IF data(ind).stat EQ 'B' THEN BEGIN
                veluga_ctree_free, data, ind, -1L, -1L, c_snap
            ENDIF ELSE BEGIN
                veluga_ctree_free, data, ind, snap_tolink, id_tolink, c_snap
			ENDELSE
        ENDIF ELSE BEGIN                ;; Tree exist
            veluga_ctree_linkbr, settings, data, ind, complete_tree, tree_key, id_tolink, snap_tolink, merit_tolink, c_snap

            IF data(ind).stat EQ 'B' THEN BEGIN
                veluga_ctree_free, data, ind, -1L, -1L, c_snap
            ENDIF ELSE BEGIN
                tkey    = snap_tolink + tree_key(0)*id_tolink
                tind2   = tree_key(kval)
                tdum    = *complete_tree(tind2)
                veluga_ctree_free, data, ind, tdum.snap(0), tdum.id(0), c_snap
            ENDELSE
        ENDELSE
    ENDFOR
END
;;-----
;; BRANCH COMPARE
;;-----
FUNCTION veluga_ctree_brcompare, settings, s0, id0, slist, idlist
	veluga 	= settings.veluga
    pid0    = veluga->r_pid(s0, id0, horg=settings.horg)
    pweight0= veluga_ctree_getweight(pid0)

    IF N_ELEMENTS(slist) EQ 1L THEN BEGIN
        pid1    = veluga->r_pid(slist(0), idlist(0), horg=settings.horg)
        pweight1= veluga_ctree_getweight(pid1)
        factor  = 1.d / settings.ctree_n_step_N
    ENDIF ELSE BEGIN
        cpid    = veluga_ctree_collectpidalongbranch(settings, slist, idlist)
        pid1    = cpid.pid
        pweight1= cpid.weight
        factor  = cpid.n_con * 1.d / settings.ctree_n_step_N
                ;pid1   = p_ctree_collectpidonbranch(tree_set, slist, idlist)
    ENDELSE
        ;; factor is to be 1 if all pid1 exist in a galaxy during n_step_bw

    larr = LONARR(20) & darr = DBLARR(20)
    ftr_name        = settings.dir_lib + '/src/fortran/get_merit.so'
        larr(0) = N_ELEMENTS(pid0)
        larr(1) = N_ELEMENTS(pid1)
        larr(2) = settings.num_thread
        larr(3) = 1L

    void    = CALL_EXTERNAL(ftr_name, 'get_merit', $
        larr, darr, pid0, pid1, pweight0, pweight1)


    RETURN, darr(0)*factor

END

;;-----
;; LINK BRANCH
;;-----
PRO veluga_ctree_linkbr, settings, data, ind, complete_tree, tree_key, idc, snapc, meritc, c_snap
	kval    = data(ind).snap + tree_key(0)*data(ind).id
    tind    = tree_key(kval)
     
    IF tind LT 0L THEN BEGIN
        veluga_ctree_makenewbr, data(ind).snap, data(ind).id, complete_tree, tree_key
        kval    = data(ind).snap + tree_key(0)*data(ind).id
        tind    = tree_key(kval)
        tmp_tree= *complete_tree(tind)
    ENDIF
    tmp_tree= *complete_tree(tind)


    kval2   = snapc + tree_key(0)*idc
    tind2   = tree_key(kval2)
    tmp_tree_toc    = *complete_tree(tind2)

    IF tmp_tree_toc.stat EQ 'MERGED' THEN BEGIN
        PRINT, 'BRANCH STOLEN'
        STOP
    ENDIF

	oldgalind       = -1L
    IF tmp_tree_toc.snap(-1) EQ tmp_tree.snap(-1) THEN BEGIN
        brorg_cut       = WHERE(tmp_tree.snap GT snapc + settings.ctree_n_step_dn, nbr_org);[0]
        brcompare_cut   = WHERE(tmp_tree_toc.snap GT snapc + settings.ctree_n_step_dn, nbr_com);[0]
        IF nbr_org EQ 0L THEN BEGIN
            brorg_cut       = WHERE(tmp_tree.snap GT snapc, nbr_org)
        ENDIF

        IF nbr_com EQ 0L THEN BEGIN
            brcompare_cut   = WHERE(tmp_tree_toc.snap GT snapc, nbr_com)
        ENDIF

        merit_com       = veluga_ctree_brcompare(settings, snapc, idc, tmp_tree_toc.snap(brcompare_cut), tmp_tree_toc.id(brcompare_cut))
        merit_org       = veluga_ctree_brcompare(settings, snapc, idc, tmp_tree.snap(brorg_cut), tmp_tree.ID(brorg_cut))

        IF merit_com GT merit_org THEN BEGIN    ;; Existing tree is better
            data(ind).stat = 'B'
            RETURN
        ENDIF
        oldgalind       = WHERE(data.id EQ tmp_tree_toc.id(0) AND data.snap EQ tmp_tree_toc.snap(0), nold)
        IF nold NE 1L THEN STOP
    ENDIF

    ;; Leave Old tree
    cut_toc = WHERE(tmp_tree_toc.snap GT snapc, nnn)
    IF nnn GE 1L THEN BEGIN
        cut_toc_merge   = WHERE(tmp_tree_toc.m_snap GT snapc, nmerge)
        IF nmerge GE 1L THEN BEGIN
            mid     = tmp_tree_toc.m_id(cut_toc_merge)
            msnap   = tmp_tree_toc.m_snap(cut_toc_merge)
            mbid    = tmp_tree_toc.m_bid(cut_toc_merge)
            mmerit  = tmp_tree_toc.m_merit(cut_toc_merge)
        ENDIF ELSE BEGIN
            mid     = -1L
            msnap   = -1L
            mbid    = -1L
            mmerit  = -1L
        ENDELSE

        IF nnn EQ 0L THEN STOP
            old_tree        = {$
                ID              : tmp_tree_toc.id(cut_toc), $
                SNAP            : tmp_tree_toc.snap(cut_toc), $
                STAT            : 'main', $
                P_SNAP          : tmp_tree_toc.p_snap(cut_toc), $
                P_ID            : tmp_tree_toc.p_id(cut_toc), $
                P_MERIT         : tmp_tree_toc.p_merit(cut_toc), $
                M_ID            : mid, $
                M_MERIT         : mmerit, $
                M_BID           : mbid, $
                M_SNAP          : msnap, $
                D_SNAP          : tmp_tree_toc.d_snap(cut_toc), $
                D_ID            : tmp_tree_toc.d_id(cut_toc), $
                ENDIND          : nnn-1L, $
                NUMPROG         : nmerge}
    ENDIF ELSE BEGIN
        old_tree        = -1L
    ENDELSE

    cut_org = WHERE(tmp_tree.snap GE data(ind).snap, norg)
    IF norg EQ 0L THEN STOP

    cut_toc = WHERE(tmp_tree_toc.snap LE snapc, ntoc)
    IF ntoc EQ 0L THEN STOP

    ;; LINK
    m_id    = [tmp_tree_toc.m_id,           tmp_tree.m_id]
    m_merit = [tmp_tree_toc.m_merit,        tmp_tree.m_merit]
    m_bid   = [tmp_tree_toc.m_bid,          tmp_tree.m_bid]
    m_snap  = [tmp_tree_toc.m_snap,         tmp_tree.m_snap]
    mcut    = WHERE(m_id GE 0L, nm)
    IF nm EQ 0L THEN BEGIN
        m_id    = [-1L]
        m_merit = [-1.0d]
        m_bid   = [-1L]
        m_snap  = [-1L]
    ENDIF ELSE BEGIN
        m_id    = m_id(mcut)
        m_merit = m_merit(mcut)
        m_bid   = m_bid(mcut)
        m_snap  = m_snap(mcut)
        ;; what if m_snap contains snapshots that are not in the tree due to pruning
    ENDELSE

    IF N_ELEMENTS(cut_org) EQ 1L THEN BEGIN
        p_snap  = [tmp_tree_toc.P_SNAP(cut_toc), snapc]
        p_id    = [tmp_tree_toc.P_ID(cut_toc), idc]
        p_merit = [tmp_tree_toc.p_merit(cut_toc), meritc]
    ENDIF ELSE BEGIN
        p_snap  = [tmp_tree_toc.P_SNAP(cut_toc), snapc, tmp_tree.p_snap(cut_org(1L:*))]
        p_id    = [tmp_tree_toc.P_ID(cut_toc), idc, tmp_tree.p_id(cut_org(1L:*))]
        p_merit = [tmp_tree_toc.p_merit(cut_toc), meritc, tmp_tree.p_merit(cut_org(1L:*))]
    ENDELSE

    IF N_ELEMENTS(cut_toc) EQ 1L THEN BEGIN
        d_snap  = [-1L, tmp_tree.snap(cut_org(0)), tmp_tree.d_snap(cut_org)]
        d_id    = [-1L, tmp_tree.id(cut_org(0)), tmp_tree.d_id(cut_org)]
    ENDIF ELSE BEGIN
        d_snap  = [tmp_tree_toc.d_snap(cut_toc(0L:-2L)), tmp_tree.snap(cut_org(0)), tmp_tree.d_snap(cut_org)]
        d_id    =  [tmp_tree_toc.d_id(cut_toc(0L:-2L)), tmp_tree.id(cut_org(0)), tmp_tree.d_id(cut_org)]
    ENDELSE

    tmp     = {$
        ID      : [tmp_tree_toc.id(cut_toc)     , tmp_tree.id(cut_org)], $
        SNAP    : [tmp_tree_toc.snap(cut_toc)   , tmp_tree.snap(cut_org)], $
        STAT    : 'main', $
        P_SNAP  : p_snap, $
        P_ID    : p_id ,$
        P_MERIT : p_merit, $
        M_ID    : m_id, $
        M_MERIT : m_merit, $
        M_BID   : m_bid, $
        M_SNAP  : m_snap, $
        D_SNAP  : d_snap, $
        D_ID    : d_id, $
        ENDIND  : N_ELEMENTS(cut_toc) + N_ELEMENTS(cut_org) - 1L, $
        NUMPROG : tmp_tree_toc.numprog + tmp_tree.numprog}

    ;; initialize
    key_vals_org    = tmp_tree.snap + tree_key(0)*tmp_tree.id
    tree_key(key_vals_org)  = -1L

    key_vals_new    = tmp_tree_toc.snap + tree_key(0)*tmp_tree_toc.id
    tree_key(key_vals_new)  = -1L

    key_vals        = tmp.snap + tree_key(0)*tmp.id
    tree_key(key_vals)      = tind
    PTR_FREE, complete_tree(tind)
    complete_tree(tind)     = PTR_NEW(tmp, /no_copy)

    IF TYPENAME(old_tree) NE 'LONG' THEN BEGIN
        key_vals_old    = old_tree.snap + old_tree.id*tree_key(0)
        tree_key(key_vals_old)  = tind2
        PTR_FREE, complete_tree(tind2)

        IF oldgalind GE 0L THEN BEGIN
            veluga_ctree_free, data, oldgalind, old_tree.snap(0), old_tree.id(0), c_snap
        ENDIF

        complete_tree(tind2)    = PTR_NEW(old_tree, /no_copy)
    ENDIF
END
;;-----
;; FREE MEMORY
;;-----
PRO veluga_ctree_free, data, ind, s_end, id_end, c_snap

    data(ind).detstat       = -1L
    PTR_FREE, data(ind).p_list
    data(ind).n_ptcl        = -1L
    IF s_end LT 0L THEN BEGIN
        data(ind).list.merit    = -1.d
        data(ind).list.id       = -1L
        data(ind).list.snap     = -1L
        data(ind).list_n        = 0L
    ENDIF ELSE IF s_end GT c_snap THEN BEGIN
        data(ind).id            = id_end
        data(ind).snap          = s_end
        cut     = WHERE(data(ind).list.snap LT s_end AND data(ind).list.snap GT 0L, ncut)
        IF ncut EQ 0L THEN BEGIN
            data(ind).list_n        = 0L
            data(ind).list.merit    = -1.d
            data(ind).list.snap     = -1L
            data(ind).list.id               = -1L
        ENDIF ELSE BEGIN
            data(ind).list(0L:ncut-1L).merit        = data(ind).list(cut).merit
            data(ind).list(0L:ncut-1L).id           = data(ind).list(cut).id
            data(ind).list(0L:ncut-1L).snap         = data(ind).list(cut).snap
            data(ind).list_n                                        = ncut

            data(ind).list(ncut:-1L).merit          = -1.d
            data(ind).list(ncut:-1L).id             = -1L
            data(ind).list(ncut:-1L).snap           = -1L
        ENDELSE
    ENDIF ELSE IF s_end LE c_snap AND s_end GT 0L THEN BEGIN
        data(ind).id            = id_end
        data(ind).snap          = s_end
        data(ind).list.merit 	= -1.d
        data(ind).list.id       = -1L
        data(ind).list.snap     = -1L
        data(ind).list_n        = 0L
    ENDIF
END

;;-----
;; MAKE NEW BRANCH
;;-----
PRO veluga_ctree_makenewbr, snap0, id0, ctree, tkey

        tmp     = {ID:[id0], SNAP:[snap0], STAT:'main', P_SNAP:[-1L], P_ID:[-1L], P_MERIT:[-1.d], $
                M_ID:[-1L], M_MERIT:[-1.d], M_BID:[-1L], M_SNAP:[-1L], $
                D_SNAP:[-1L], D_ID:[-1L], ENDIND:0L, NUMPROG:0L}
        tmp_ptr = PTR_NEW(tmp, /no_copy)
        ctree   = [ctree, tmp_ptr]
        tind    = N_ELEMENTS(ctree)-1L
        tkey(snap0 + tkey(0)*id0)       = tind
END

;;-----
;; EXPAND BRANCH
;;-----
PRO veluga_ctree_expandbr, data, ind, complete_tree, tree_key, idc, snapc, meritc
    kval    = data(ind).snap + tree_key(0)*data(ind).id
    tind    = tree_key(kval)
    
    IF tind LT 0L THEN BEGIN
        veluga_ctree_makenewbr, data(ind).snap, data(ind).id, complete_tree, tree_key
        tind    = tree_key(data(ind).snap + tree_key(0)*data(ind).id)
    ENDIF

    tmp_tree= *complete_tree(tind)

    IF N_ELEMENTS(tmp_tree.id) EQ 1L THEN BEGIN
        p_snap  = [-1L, snapc]
        p_id    = [-1L, idc]
        p_merit = [-1.d, meritc]
    ENDIF ELSE BEGIN
        p_snap  = [-1L, snapc, tmp_tree.p_snap(1L:*)]
        p_id    = [-1L, idc , tmp_tree.p_id(1L:*)]
        p_merit = [-1.d, meritc, tmp_tree.p_merit(1L:*)]
    ENDELSE

    tmp     = {$
        ID      : [idc, tmp_tree.id], $
        SNAP    : [snapc, tmp_tree.snap], $
        STAT    : 'main', $
        P_SNAP  : p_snap, $
        P_ID    : p_id, $
        P_MERIT : p_merit, $
        M_ID    : tmp_tree.m_id, $
        M_MERit : tmp_tree.m_merit, $
        M_BID   : tmp_tree.m_bid, $
        M_SNAP  : tmp_tree.m_snap, $
        D_SNAP  : [tmp_tree.snap(0), tmp_tree.d_snap], $
        D_ID    : [tmp_tree.id(0), tmp_tree.d_id], $
        ENDIND  : tmp_tree.endind + 1L, $
        NUMPROG : tmp_tree.numprog}

    kval_new        = snapc + tree_key(0)*idc
    
    IF tree_key(kval_new) GE 0L THEN BEGIN
        PRINT, 'tree corrupted'
        STOP
    ENDIF

    tree_key(kval_new)      = tind
    PTR_FREE, complete_tree(tind)
    complete_tree(tind)     = PTR_NEW(tmp, /no_copy)
END

;;-----
;; MAIN
;;	Complete tree corrector made by Jinsu Rhee
;;-----
PRO veluga_ctree, header, num_thread=num_thread, horg=horg

	
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
	settings   = CREATE_STRUCT(settings, 'dir_tree', dir_tree)

	TIC & TOC

	snapdata 	= veluga_ctree_getsnapinfo(settings, veluga)


	settings	= CREATE_STRUCT(settings, 'veluga', veluga, 'sinfo', snapdata.sinfo, 'slist', snapdata.slist)
	
	;;-----
	;; LOAD the original tree
	;;-----
	RESTORE, settings.dir_tree + '/tfout/tree.sav'

	;;-----
	;; ALLOCATE
	;;-----
	gal 	= veluga->r_gal(MAX(settings.ctree_snap), -1L, horg=horg, GProp=['ID'])
	data 	= veluga_ctree_allocate(settings, complete_tree, tree_key, N_ELEMENTS(gal))

	ng0 	= 0L
	ng1 	= N_ELEMENTS(gal)-1L

	veluga_ctree_intputgal, settings, complete_tree, tree_key, data, gal, ng0, ng1



	;;-----
	;; MAIN LOOP
	;;		TO DO
	;;			Merit calculation considering weights
	;;-----
	FOR i=N_ELEMENTS(settings.slist)-1L, 0L, -1L DO BEGIN
		
		c_snap 	= settings.slist(i)

		PRINT, '%123123-----'
        PRINT, ''
        PRINT, '                TREE CONNECTION AT SNAP = ' + STRING(c_snap,format='(I4.4)') + ' ( ' + STRING(N_ELEMENTS(data),format='(I6)') + ' gals )'

      

        ;;-----
        ;; CLASSIFY GALAXIES WITH BROKEN TREES
        ;;-----
        TIC
        veluga_ctree_classify, settings, data, c_snap, number
        TOC, elapsed_time=t_classify
        IF number.T EQ N_ELEMENTS(data) THEN BEGIN
        	PRINT, '                        SKIP due to all galaxies having trees'
        	CONTINUE
        ENDIF

        ;;-----
        ;; Determine End point
        ;;-----
        TIC
        velgua_ctree_detend, settings, data, complete_tree, tree_key
       	TOC, elapsed_time=t_detend

       	;;-----
       	;; Collect particles of galaxies for their merit to be computed
       	;;-----
       	TIC
       	pid 	= veluga_ctree_collectpid(settings, data, complete_tree, tree_key)
       	TOC, elapsed_time=t_cpid

       	;;-----
        ;; Read particles at this snapshot
        ;;-----
        TIC
        pid0    = veluga_ctree_readsnap(settings, data, c_snap)
        TOC, elapsed_time=t_rsnap

        ;;-----
        ;; Merit Calcultion
        ;;-----
        TIC
        veluga_ctree_commerit, settings, data, pid, pid0, c_snap
        TOC, elapsed_time=t_merit

        ;;-----
        ;; LINK TREE
        ;;-----
        TIC
        veluga_ctree_link, settings, data, number, c_snap, complete_tree, tree_key
        TOC, elapsed_time=t_link


        ;;-----
        ;; TODO
        ;;-----

        ; If tree is connected or renewed, set detstat to be -1
        ;               free p_list and set n_ptcl to be negative
        ;               free list and list_n

        PRINT, '                Time report [sec]'
        PRINT, '                        Classify Galaxies :', t_classify
        PRINT, '                        Det- Branch End   :', t_detend
        PRINT, '                        Collect PID       :', t_cpid
        PRINT, '                        Read Snap ptcls   :', t_rsnap
        PRINT, '                        Compute Merits    :', t_merit
        PRINT, '                        Link Branch       :', t_link

        IF i MOD 5L EQ 0L THEN $
            SAVE, FILENAME=settings.dir_tree + '/tfout/ctree_' + STRING(c_snap,format='(I4.4)') + '.sav', settings, data, c_snap, complete_tree, tree_key

        IF c_snap LE settings.ctree_snap(0) THEN BEGIN
            veluga_ctree_classify, settings, data, c_snap, number
            veluga_ctree_detend, settings, data, complete_tree, tree_key
            REPEAT BEGIN
                veluga_ctree_link, settings, data, number, c_snap, complete_tree, tree_key
                veluga_ctree_classify, settings, data, c_snap, number
                veluga_ctree_detend, settings, data, complete_tree, tree_key

            ENDREP UNTIL MAX(data.list_n) EQ 0L
            BREAK
        ENDIF
    ENDFOR

    SAVE, filename=settings.dir_tree + '/tfout/ctree.sav', complete_tree, tree_key
    SAVE, filename=settings.dir_tree + '/tfout/ctree_dat.sav', settings, data
    PRINT, 'Done ^-^'

END