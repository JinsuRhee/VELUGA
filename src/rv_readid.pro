FUNCTION rv_readid, settings, veluga, runstat, run=run

dir_data        = runstat.dir
;;-----
;; Check procedure set
;;-----
IF run EQ 0L THEN RETURN, PTR_NEW({p_id:-1, b_ind:-1, u_ind:-1},/no_copy)
IF run EQ 1L THEN BEGIN
	IF STRLEN(FILE_SEARCH(dir_data + 'rv_id.sav')) GE 5L THEN BEGIN
		RESTORE, dir_data + 'rv_id.sav'
		RETURN, PTR_NEW(output,/no_copy)
	ENDIF ELSE BEGIN
		run	= 2L
	ENDELSE
ENDIF
IF run EQ 2L THEN BEGIN
	;PRINT, '        %%%%% (No previous works are found)'

	;;-----
	;; I/O Settings
	;;-----
	IF settings.horg EQ 'h' THEN pref = '/halo'
	IF settings.horg EQ 'g' THEN pref = '/galaxy'

	dum_fname_pt	= dir_data + STRTRIM(pref,2) + '.dat.catalog_groups.*'
	dum_fname_ptb	= dir_data + STRTRIM(pref,2) + '.dat.catalog_particles.*'
	dum_fname_ptu	= dir_data + STRTRIM(pref,2) + '.dat.catalog_particles.unbound.*'

	dum_fname	= FILE_SEARCH(dum_fname_pt)
	dum_fname_bdn	= FILE_SEARCH(dum_fname_ptb)
	dum_fname_ubd	= FILE_SEARCH(dum_fname_ptu)

	;;-----
	;; Read Ptcl ID Info
	;;-----

        str     = ' '
        OPENR, 10, dum_fname(0)
        READF, 10, str
        str2    = STRSPLIT(str, ' ', /extract)
        n_mpi	= LONG(str2(1))
        CLOSE, 10

        n_obj   = LONARR(n_mpi)
        FOR i=0L, n_mpi-1L DO BEGIN
                OPENR, 10, dum_fname(i)
                READF, 10, str & READF, 10, str
                str2    = STRSPLIT(str, ' ', /extract)
                n_obj(i)= LONG(str2(0))
                IF i EQ 0L THEN n_tot = LONG(str2(1))
                CLOSE, 10
        ENDFOR

        n_part  = LONARR(n_tot,4)
        i0      = 0L
        FOR i=0L, n_mpi-1L DO BEGIN
                IF n_obj(i) EQ 0L THEN CONTINUE
                OPENR, 10, dum_fname(i)
                READF, 10, str & READF, 10, str

                n_part(i0:i0 + n_obj(i)-1, 0) = i               ;; NUM MPI
                FOR j=0L, n_obj(i)-1L DO BEGIN
                        READF, 10, str
                        n_part(i0 + j,1) = LONG(str)            ;; # of ALL PTCLS
                ENDFOR
                FOR j=0L, n_obj(i)-1L DO BEGIN
                        READF, 10, str
                        n_part(i0 + j,2) = LONG(str)            ;; BDN INDICES
                ENDFOR
                FOR j=0L, n_obj(i)-1L DO BEGIN
                        READF, 10, str
                        n_part(i0 + j,3) = LONG(str)            ;; UBD INDICES
                ENDFOR

                i0 = i0 + n_obj(i)
                CLOSE, 10
        ENDFOR

	i0 = 0L & j0 = 0L & k0 = 0L
	l0 = 0L & m0 = 0L

        FOR i=0L, n_mpi-1L DO BEGIN
                OPENR, 10, dum_fname_bdn(i)
                OPENR, 11, dum_fname_ubd(i)
                READF, 10, str & READF, 10, str
                IF i EQ 0L THEN BEGIN
                        str2    = STRSPLIT(str, ' ', /extract)
                        n_bdn_tot= LONG(str2(1))
                ENDIF
                READF, 11, str & READF, 11, str
                IF i EQ 0L THEN BEGIN
                        str2    = STRSPLIT(str, ' ', /extract)
                        n_ubd_tot= LONG(str2(1))
                ENDIF

                IF i EQ 0L THEN BEGIN
                        id_bdn  = LON64ARR(n_bdn_tot + 1L)
                        id_ubd  = LON64ARR(n_ubd_tot + 1L)
                        bdn_ind = LONARR(n_tot,2)
                        ubd_ind = LONARR(n_tot,2)
                ENDIF

                FOR j=0L, n_obj(i)-1L DO BEGIN
                        IF j NE n_obj(i) -1L THEN BEGIN
                                n_tot = n_part(i0 + j,1)
                                n_bnd = n_part(i0 + j + 1,2) - n_part(i0 + j,2)
                                n_ubd = n_part(i0 + j + 1,3) - n_part(i0 + j,3)
                        ENDIF ELSE BEGIN
                                n_tot = n_part(i0 + j,1)
                                n_bnd = FILE_LINES(dum_fname_bdn(i)) - 2L - n_part(i0 + j,2)
                                n_ubd = FILE_LINES(dum_fname_ubd(i)) - 2L - n_part(i0 + j,3)
                        ENDELSE

                        IF n_bnd GE 1L THEN BEGIN
                                bdn_ind(i0 + j,0) = j0
				id_bdn_dum	= LON64ARR(n_bnd)
                                FOR k=0L, n_bnd-1L DO BEGIN
                                        READF, 10, str
					id_bdn_dum(k) = LONG64(str)
                                ENDFOR
				;js_makearr, id_bdn, id_bdn_dum, j0, unitsize=100000L, type='L64'
                                veluga->g_makearr, id_bdn, id_bdn_dum, j0, unitsize=100000L, type='L64'
                                bdn_ind(i0 + j,1) = j0 - 1L
                        ENDIF ELSE BEGIN
                                bdn_ind(i0 + j,0) = j0
                                bdn_ind(i0 + j,1) = j0

				IF j0 EQ 0L THEN id_bdn(0) = -1L
				IF j0 GE 1L THEN veluga->g_makearr, id_bdn, [-9223372036854775808], $
					j0, unitsize=100000L, type='L64'
                        ENDELSE

                        IF n_ubd GE 1L THEN BEGIN
                                ubd_ind(i0 + j,0) = k0
				id_ubd_dum	= LON64ARR(n_ubd)
                                FOR k=0L, n_ubd-1L DO BEGIN
                                        READF, 11, str
					id_ubd_dum(k) = LONG64(str)
                                ENDFOR
				;js_makearr, id_ubd, id_ubd_dum, k0, unitsize=100000L, type='L64'
                                veluga->g_makearr, id_ubd, id_ubd_dum, k0, unitsize=100000L, type='L64'
                                ubd_ind(i0 + j,1) = k0 - 1L
                        ENDIF ELSE BEGIN
                                ubd_ind(i0 + j,0) = k0
                                ubd_ind(i0 + j,1) = k0
				IF k0 EQ 0L THEN id_ubd(0) = -1L
				IF k0 GE 1L THEN veluga->g_makearr, id_ubd, [-9223372036854775808], $
					k0, unitsize=100000L, type='L64'	;;;
                        ENDELSE
                ENDFOR
  
                CLOSE, 10 & CLOSE, 11
                i0 = i0 + n_obj(i)
        ENDFOR
 
        id_bdn = id_bdn(0L:j0-1L) & id_ubd = id_ubd(0L:k0-1L)

	output  = CREATE_STRUCT('p_id', [id_bdn, id_ubd])
	output  = CREATE_STRUCT(output, 'b_ind', bdn_ind, 'u_ind', ubd_ind + n_elements(id_bdn))

	IF settings.pp_saveprocess EQ 1L THEN SAVE, filename=dir_data + 'rv_id.sav', output
	RETURN, PTR_NEW(output,/no_copy)
ENDIF
END
