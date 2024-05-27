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


PRO veluga_makebr, header, num_thread=num_thread, horg=horg

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

	treeset 	= veluga_makebr_init(settings, veluga)

	;;-----
	;; ALLOCATE
	;;-----
	
	;runstat 	= REPLICATE({snap:-1L, iscatalog:-1L, dir:' ', elt:DBLARR(6), $
	;	rv_raw:PTR_NEW(1), rv_id:PTR_NEW(1), rv_ptmatch:PTR_NEW(1), rv_bprop:PTR_NEW(1) $
	;		}, settings.pp_snap(1)-settings.pp_snap(0)+1L)


	;;-----
	;; COMPILE ALL PROCEDURES FIRST
	;;-----
	;veluga__initcompile, settings, veluga, runstat(0)

	;;-----
	;;
	;;-----
END