PRO veluga_pp_runcheck, settings, runstat, ind

	dir 	= settings.dir_catalog
	IF settings.horg EQ 'g' THEN dir += 'Galaxy/'
	IF settings.horg EQ 'h' THEN dir += 'Halo/'
	dir 	+= 'snap_' + STRING(runstat(ind).snap, format='(I4.4)')
	runstat(ind).dir 	= dir
	fname 	= dir + '/*.dat.properties.0'
	isfile 	= STRLEN(FILE_SEARCH(fname))
	IF isfile GE 5L THEN runstat(ind).iscatalog = 1L
	RETURN
END

PRO veluga_pp, header, num_thread=num_thread, horg=horg

	;;-----
	;; INITIAL SETTINGS
	;;-----
	IF ~KEYWORD_SET(num_thread) THEN num_thread = 10L
	IF ~KEYWORD_SET(horg) THEN horg = 'g'

	settings 	= CREATE_STRUCT(settings, 'horg', horg)
	;;-----
	;; CALL OBJECT
	;;-----
	veluga	= OBJ_NEW('veluga', header, num_thread=num_thread)

	settings 	= veluga->getheader()

	;;-----
	;; COMPILE ALL PROCEDURES FIRST
	;;-----
	;void 	= rv_RawCatalog(settings, ' ', run=0L)
	;void 	= rv_ReadID(settings, ' ', run=0L)
	;void 	= rv_PTMatch(settings, ' ', run=0L)
	;void 	= rv_GProp(settings, ' ', run=0L)
	;void 	= rv_RawCatalog(settings, ' ', run=0L)
	;void 	= rv_RawCatalog(settings, ' ', run=0L)
	;void 	= rv_RawCatalog(settings, ' ', run=0L)
	;void 	= r

	;;-----
	;; RUN
	;;-----
	runstat 	= REPLICATE({snap:-1L, iscatalog:-1L, dir:' ', $
		rv_raw:PTR_NEW(1), rv_id:PTR_NEW(1), rv_ptmatch:PTR_NEW(1), rv_bprop:PTR_NEW(1) $
			}, settings.pp_snap(1)-settings.pp_snap(0)+1L)

	ind 		= 0L
	FOR i=settings.pp_snap(0), settings.pp_snap(1), settings.pp_snap(2) DO BEGIN

		
		TIC
		runstat(ind).snap 	= i
		veluga->ppout, '		snapshot ' + STRING(i,format='(I4.4)') + ' is starting'

		;;----- RUN CHECK
		veluga_pp_runcheck, settings, runstat, ind

		;;-----
		;; READ Raw catalog data
		;;-----
		veluga->ppout2, 'Reading The Raw Catalog...'
		runstat(ind).rv_raw 	= rv_RawCatalog(settings, veluga, runstat(ind), run=settings.pp_runtype.catalog)

		;;-----
		;; READ Particle IDs
		;;-----
		veluga->ppout2, 'Reading The member ptcl ID...'
		runstat(ind).rv_id 	= rv_ReadID(settings, veluga, runstat(ind), run=settings.pp_runtype.ptclid)

		;;-----
		;; PTCL Matching
		;;-----
		veluga->ppout2, 'Particle membership matching'
		runstat(ind).rv_ptmatch 	= rv_PTMatch(settings, veluga, runstat(ind), run=settings.pp_runtype.ptmatch)
		
		;;-----
		;; Bulk Properties
		;;-----
		veluga->ppout2, 'Bulk properties computations'
		runstat(ind).rv_bprop		= rv_Bprop(settings, veluga, runstat(ind), run=settings.pp_runtype.bprop)
		
		;;---- RUN STAT CHECK
		;;123123 no catalog file
		;;		make txt file maybe?
		TOC, elapsed_time = elt
		ind ++
	ENDFOR

END
