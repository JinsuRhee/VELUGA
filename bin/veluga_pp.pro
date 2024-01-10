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
PRO veluga_pp_initcompile, settings, veluga, runstat
	void 	= rv_RawCatalog(settings, veluga, runstat(0), run=0L)
	void 	= rv_ReadID(settings, veluga, runstat(0), run=0L)
	void 	= rv_PTMatch(settings, veluga, runstat(0), run=0L)
	void 	= rv_BProp(settings, veluga, runstat(0), run=0L)
	rv_save, settings, veluga, runstat(0), run=0L

	simple_write_hdf5, 1, 1, 1, /skip
	dumarr 	= FINDGEN(10,10)
	dumind 	= ARRAY_INDICES(dumarr, [1])
	TIC
	TOC, elapsed_time=elt
END

PRO veluga_pp, header, num_thread=num_thread, horg=horg

	;;-----
	;; INITIAL SETTINGS
	;;-----
	IF ~KEYWORD_SET(num_thread) THEN num_thread = 10L
	IF ~KEYWORD_SET(horg) THEN horg = 'g'


	;;-----
	;; CALL OBJECT
	;;-----
	veluga	= OBJ_NEW('veluga', header, num_thread=num_thread)

	settings 	= veluga->getheader()
	settings 	= CREATE_STRUCT(settings, 'horg', horg)

	runstat 	= REPLICATE({snap:-1L, iscatalog:-1L, dir:' ', elt:DBLARR(6), $
		rv_raw:PTR_NEW(1), rv_id:PTR_NEW(1), rv_ptmatch:PTR_NEW(1), rv_bprop:PTR_NEW(1) $
			}, settings.pp_snap(1)-settings.pp_snap(0)+1L)


	;;-----
	;; COMPILE ALL PROCEDURES FIRST
	;;-----
	veluga_pp_initcompile, settings, veluga, runstat(0)


	;void 	= rv_RawCatalog(settings, ' ', run=0L)
	;void 	= rv_RawCatalog(settings, ' ', run=0L)
	;void 	= rv_RawCatalog(settings, ' ', run=0L)
	;void 	= r

	;;-----
	;; RUN
	;;-----
	

	ind 		= 0L
	FOR i=settings.pp_snap(0), settings.pp_snap(1), settings.pp_snap(2) DO BEGIN

		
		
		runstat(ind).snap 	= i
		veluga->ppout, '		snapshot ' + STRING(i,format='(I4.4)') + ' is starting'

		;;----- RUN CHECK
		veluga_pp_runcheck, settings, runstat, ind

		;;-----
		;; READ Raw catalog data
		;;-----
		TIC
		veluga->ppout2, 'Reading The Raw Catalog...'
		runstat(ind).rv_raw 	= rv_RawCatalog(settings, veluga, runstat(ind), run=settings.pp_runtype.load_catalog)
		TOC, elapsed_time=elt_lc
		veluga->ppout2, '(Done in' + STRING(elt_lc,format='(F9.4)') + ' [s])'
		veluga->ppout2, ' '
		runstat(ind).elt(1) = elt_lc

		;;-----
		;; READ Particle IDs
		;;-----
		veluga->ppout2, 'Reading The member ptcl ID...'
		TIC
		runstat(ind).rv_id 	= rv_ReadID(settings, veluga, runstat(ind), run=settings.pp_runtype.read_ptclid)
		TOC, elapsed_time=elt_rp
		veluga->ppout2, '(Done in' + STRING(elt_rp,format='(F9.4)') + ' [s])'
		veluga->ppout2, ' '
		runstat(ind).elt(2) = elt_rp

		;;-----
		;; PTCL Matching
		;;-----
		veluga->ppout2, 'Particle membership matching'
		TIC
		runstat(ind).rv_ptmatch 	= rv_PTMatch(settings, veluga, runstat(ind), run=settings.pp_runtype.member_match)
		TOC, elapsed_time=elt_mm
		veluga->ppout2, '(Done in' + STRING(elt_mm,format='(F9.4)') + ' [s])'
		veluga->ppout2, ' '
		runstat(ind).elt(3) = elt_mm

		;;-----
		;; Bulk Properties
		;;-----
		veluga->ppout2, 'Bulk properties computations'
		TIC
		runstat(ind).rv_bprop		= rv_Bprop(settings, veluga, runstat(ind), run=settings.pp_runtype.compute_bulk)
		TOC, elapsed_time=elt_cb
		veluga->ppout2, '(Done in' + STRING(elt_cb,format='(F9.4)') + ' [s])'
		veluga->ppout2, ' '
		runstat(ind).elt(4) = elt_cb

		;;-----
		;; SAVE catalog
		;;-----
		veluga->ppout2, 'Making hdf5 catalog'
		TIC
		rv_save, settings, veluga, runstat(ind), run=settings.pp_runtype.save
		TOC, elapsed_time=elt_sv
		veluga->ppout2, '(Done in' + STRING(elt_sv,format='(F9.4)') + ' [s])'
		veluga->ppout2, ' '
		runstat(ind).elt(5) = elt_sv

		;;-----
		;; FREE MEMORY
		;;-----
		PTR_FREE, runstat(ind).rv_raw
		PTR_FREE, runstat(ind).rv_id
		PTR_FREE, runstat(ind).rv_ptmatch
		PTR_FREE, runstat(ind).rv_bprop

		runstat(ind).elt(0)	= TOTAL(runstat(ind).elt(1L:-1L))
		;;-----
		;; REPORT
		;;-----
		veluga->ppout2, '		snapshot ' + STRING(i,format='(I4.4)') + ' is done'
		veluga->ppout2, ' '
		veluga->ppout2, '	Total Wall clock time : ' + STRING(runstat(ind).elt(0), format='(F9.4)')
		veluga->ppout2, ' '
		
		ind ++
	ENDFOR

END
