PRO rv_save, settings, veluga, runstat, run=run

dir_data	= runstat.dir

;;-----
;; Check procedure set
;;-----
IF run EQ 0L THEN RETURN

	;;-----
	;; Make a Directory?
	;;----
	IF settings.horg EQ 'g' THEN dum = 'Galaxy/VR_Galaxy'
	IF settings.horg EQ 'h' THEN dum = 'Halo/VR_Halo'

	isfile	= STRLEN(FILE_SEARCH(settings.dir_catalog+dum))
	IF isfile LE 5L THEN SPAWN, 'mkdir ' + settings.dir_catalog + dum


	n_snap 	= runstat.snap
	fname	= settings.dir_catalog + dum + '/snap_' + STRING(n_snap,format='(I4.4)') + '.hdf5'

	;;-----
	;; Create HDF
	;;-----

	;;	----- PTCLs
	;;		POS: in Kpc
	;;		VEL: in km/s
	;;		AGE: in conformal
	;;		SF: in scale factor
	;;		GYR: in Gyr
	;;		Mass: in Solar mass

	ngal	= N_ELEMENTS((*runstat.rv_raw).id)
	nsfr 	= N_ELEMENTS(settings.SFR_R)
	nmpi 	= settings.ndomain
	nconf 	= N_ELEMENTS(settings.CONF_R)

	;;-----
	;; SAVE HDF5 by c
	;;-----

	;;----- merge catalog data
	cat_dtype	= LONARR(N_ELEMENTS(settings.column_list))

	int_dtypelist	= ['ID', 'ID_mbp', 'hostHaloID', 'numSubStruct', 'Structuretype', 'npart']
	FOR i=0L, N_ELEMENTS(settings.column_list)-1L DO BEGIN
		tmp	= settings.column_list(i)

		cut 	= WHERE(int_dtypelist EQ tmp, nislong)
		IF nislong GE 1L THEN BEGIN
			cat_dtype(i) = 1L
		ENDIF ELSE BEGIN
			cat_dtype(i) = -1L
		ENDELSE
	ENDFOR

	clist_long	= settings.column_list(WHERE(cat_dtype EQ 1L))
	clist_dbl	= settings.column_list(WHERE(cat_dtype EQ -1L))

	cat_double	= DBLARR(ngal * N_ELEMENTS(clist_dbl))
	cat_long	= LONARR(ngal * N_ELEMENTS(clist_long))

	i0	= 0L
	FOR i=0L, N_ELEMENTS(clist_long)-1L DO BEGIN
		i1	= i0 + ngal-1L

		str	= 'cat_long(' + STRTRIM(i0) + ':' + STRTRIM(i1) + ') = (*runstat.rv_raw).' + clist_long(i)
		void	= EXECUTE(str)
		i0	= i1 + 1L
	ENDFOR

	i0	= 0L
	FOR i=0L, N_ELEMENTS(clist_dbl)-1L DO BEGIN
		i1	= i0 + ngal-1L

		str	= 'cat_double(' + STRTRIM(i0) + ':' + STRTRIM(i1) + ') = (*runstat.rv_raw).' + clist_dbl(i)
		void	= EXECUTE(str)
		i0	= i1 + 1L
	ENDFOR

	;;----- Bulk properties
	bprop	= *runstat.rv_bprop
	IF N_ELEMENTS(bprop.sfr) GE 2L THEN BEGIN
		b_sfr	= TRANSPOSE(bprop.sfr)
	ENDIF ELSE BEGIN
		b_sfr	= [-1.d]
	ENDELSE

	IF N_ELEMENTS(bprop.abmag) GE 2L THEN BEGIN
		b_mag	= DBLARR(ngal*N_ELEMENTS(settings.flux_list)*N_ELEMENTS(settings.mag_r))
		i0	= 0L
		FOR i=0L, ngal-1L DO BEGIN
			i1	= i0 + N_ELEMENTS(settings.flux_list)*N_ELEMENTS(settings.mag_r)-1L
			tmp = bprop.abmag(i)
			str	= 'dummy =['
			FOR j=0L, N_ELEMENTS(settings.flux_list)-1L DO BEGIN
				str += 'tmp.' + STRTRIM(settings.flux_list(j))
				IF j LT N_ELEMENTS(settings.flux_list)-1L THEN str += ', '
			ENDFOR
			str	+= ']'
			void	= EXECUTE(str)

			b_mag(i0:i1)	= dummy
			i0	= i1 + 1L
		ENDFOR
	ENDIF ELSE BEGIN
		b_mag = [-1.d]
	ENDELSE

	IF N_ELEMENTS(bprop.sb) GE 2L THEN BEGIN
		b_sb	= DBLARR(ngal*N_ELEMENTS(settings.flux_list)*N_ELEMENTS(settings.mag_r))
		i0	= 0L
		FOR i=0L, ngal-1L DO BEGIN
			i1	= i0 + N_ELEMENTS(settings.flux_list)*N_ELEMENTS(settings.mag_r)-1L
			tmp = bprop.sb(i)
			str	= 'dummy =['
			FOR j=0L, N_ELEMENTS(settings.flux_list)-1L DO BEGIN
				str += 'tmp.' + STRTRIM(settings.flux_list(j))
				IF j LT N_ELEMENTS(settings.flux_list)-1L THEN str += ', '
			ENDFOR
			str	+= ']'
			void	= EXECUTE(str)

			b_sb(i0:i1)	= dummy
			i0	= i1 + 1L
		ENDFOR
	ENDIF ELSE BEGIN
		b_sb = [-1.d]
	ENDELSE

	b_cm	= bprop.confrac_m.aper
	b_cn	= bprop.confrac_n.aper

	IF N_ELEMENTS(bprop.isclump) GE 2L THEN BEGIN
		b_isclump	= bprop.isclump
	ENDIF ELSE BEGIN
		b_isclump 	= [-1L]
	ENDELSE

	pprop	= (*runstat.rv_ptmatch)
	IF N_ELEMENTS(pprop.dom_list) GE 2L THEN BEGIN
		b_domlist	= TRANSPOSE(pprop.dom_list)
	ENDIF ELSE BEGIN
		b_domlist	= [-1L]
	ENDELSE

	;;----- Send TO C
	ftr_name        = settings.dir_lib + 'src/fortran/save_cat_conly.so'
                larr = LONARR(100) & darr = DBLARR(20)
                larr(0) = ngal
                larr(1) = nmpi

		larr(10)	= N_ELEMENTS(settings.sfr_r)
		larr(11)	= N_ELEMENTS(settings.mag_r)
		larr(12)	= N_ELEMENTS(settings.conf_r)
		larr(13)	= N_ELEMENTS(settings.flux_list)
		larr(14)	= N_ELEMENTS((*runstat.rv_id).b_ind)
		larr(15)	= N_ELEMENTS(clist_long)
		larr(16)	= N_ELEMENTS(clist_dbl)
		larr(17)	= N_ELEMENTS(b_sfr)
		larr(18)	= N_ELEMENTS(b_mag)
		larr(19)	= N_ELEMENTS(b_sb)
		larr(20)	= N_ELEMENTS(b_isclump)
		larr(21)	= N_ELEMENTS(b_domlist)

		darr(0)		= (*runstat.rv_ptmatch).a_exp
	
        void    = CALL_EXTERNAL(ftr_name, 'save_cat', $
                larr, darr, fname, $
		DOUBLE(settings.sfr_r), DOUBLE(settings.sfr_t), DOUBLE(settings.mag_r), DOUBLE(settings.conf_r), settings.flux_list, $ ;; send info (3-7)
		DOUBLE((*runstat.rv_raw).mass_tot),  $ 		;; bulk 8
		DOUBLE((*runstat.rv_raw).r_halfmass),  $ 	;; bulk 9
		DOUBLE((*runstat.rv_raw).mvir),  $ 		;; bulk 10 
		DOUBLE((*runstat.rv_raw).rvir),  $ 		;; bulk 11
		DOUBLE((*runstat.rv_raw).mass_200crit),  $	;; bulk 12
		DOUBLE((*runstat.rv_raw).r_200crit),  $		;; bulk 13
		LONG((*runstat.rv_raw).ID),  $			;; bulk 14
		DOUBLE(TRANSPOSE((*runstat.rv_bprop).confrac_m.aper)),  $		;; bulk 15
		DOUBLE(TRANSPOSE((*runstat.rv_bprop).confrac_n.aper)),  $		;; bulk 16
		LONG64((*runstat.rv_id).b_ind), $		;; bulk 17
		LONG64((*runstat.rv_id).u_ind), $		;; bulk 18
		LONG64((*runstat.rv_id).p_id), $		;; bulk 19
		LONG(cat_long), $				;; bulk 20
		clist_long, $					;; bulk 21
		DOUBLE(cat_double), $				;; bulk 22
		clist_dbl, $					;; bulk 23
		DOUBLE(b_sfr), $				;; bulk 24
		DOUBLE(b_mag), $				;; bulk 25
		DOUBLE(b_sb), $					;; bulk 26
		DOUBLE(b_cm), $					;; bulk 27
		DOUBLE(b_cn), $					;; bulk 28
		LONG(b_isclump), $				;; bulk 29
		LONG(b_domlist), $				;; bulk 30
		1L $ ; for dummy
		)

	;;CONF_M, CONF_N

	RETURN

	;;----- Script below is an old one
	;;-----
	;; Open HDF5
	;;-----

	fid	= h5f_create(fname)

	;;----- Write General Information
	simple_write_hdf5, settings.flux_list, 'Flux_List', 	fid

	simple_write_hdf5, settings.SFR_R, 'SFR_R',		fid
	simple_write_hdf5, settings.SFR_T, 'SFR_T', 		fid
	simple_write_hdf5, settings.MAG_R, 'MAG_R', 		fid
	simple_write_hdf5, settings.CONF_R, 'CONF_R', 			fid

	;;----- Write some bulk properties
	simple_write_hdf5, (*runstat.rv_raw).mass_tot, 'Mass_tot', fid
	simple_write_hdf5, (*runstat.rv_raw).r_halfmass, 'R_HalfMass', fid
	simple_write_hdf5, (*runstat.rv_raw).mvir, 'Mvir', fid
	simple_write_hdf5, (*runstat.rv_raw).rvir, 'Rvir', fid
	simple_write_hdf5, (*runstat.rv_raw).id, 'ID', fid
	simple_write_hdf5, TRANSPOSE((*runstat.rv_bprop).confrac_m.aper), 'CONF_M', fid
	simple_write_hdf5, TRANSPOSE((*runstat.rv_bprop).confrac_n.aper), 'CONF_N', fid

	FOR i=0L, ngal - 1L DO BEGIN

		ib = -1L & iu = -1L & ptcl_id = -1L
		IF N_ELEMENTS((*runstat.rv_id).b_ind) GE 2 THEN BEGIN
			ib = (*runstat.rv_id).b_ind(i,*) & iu = (*runstat.rv_id).u_ind(i,*)
			ptcl_id	= [(*runstat.rv_id).p_id(ib(0):ib(1)), $
				(*runstat.rv_id).p_id(iu(0):iu(1))]

			cut	= WHERE(ptcl_id GT -922337203685477580LL, ncut)
			IF ncut NE (*runstat.rv_raw).npart(i) THEN BEGIN
				PRINT, 'WRONG PTCL ID ARRAY'
				STOP
			ENDIF
			ptcl_id	= ptcl_id(cut)
		ENDIF

		;;----- Create Groups for this galaxy
		idstr	= 'ID_' + STRING((*runstat.rv_raw).id(i), format='(I6.6)')
		gpstr	= idstr + '/G_Prop'
		ppstr	= idstr + '/P_Prop'

		void	= H5G_CREATE(fid, idstr) 
		void	= H5G_CREATE(fid, gpstr)
		void	= H5G_CREATE(fid, ppstr)


		;----- Write Raw Catalog Properties
		FOR j=0L, N_ELEMENTS(settings.column_list)-1L DO BEGIN
			str	= 'tmp = [(*runstat.rv_raw).' + settings.column_list(j) + '(i)]'
			void	= EXECUTE(str)
			simple_write_hdf5, tmp, + gpstr + '/G_' + settings.column_list(j),	fid
		ENDFOR

		;----- Write Bulk Properties
		bprop 	= *runstat.rv_bprop

		sfr = -1.d
		IF N_ELEMENTS(bprop.sfr) GE 2L THEN sfr = REFORM(bprop.sfr(i,*), nsfr)
		simple_write_hdf5, sfr, gpstr + '/G_SFR', fid

		
		FOR fi=0L, N_ELEMENTS(settings.flux_list)-1L DO BEGIN
			nullarr	= -1.d
			txtdum	= '/G_ABmag_' + STRTRIM(settings.flux_list(fi),2)

			IF N_ELEMENTS(bprop.abmag) GE 2L THEN void2	= EXECUTE('nullarr = bprop.abmag(' + STRTRIM(i,2) + ').' + settings.flux_list(fi))
			
			simple_write_hdf5, nullarr, gpstr + txtdum, fid
		ENDFOR
		
		;nuv = -1.d & mu = -1.d & mg = -1.d & mr = -1.d & mi = -1.d & mz = -1.d
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN nuv = bprop.abmag(i).nuv
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mu = bprop.abmag(i).u
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mg = bprop.abmag(i).g
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mr = bprop.abmag(i).r
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mi = bprop.abmag(i).i
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mz = bprop.abmag(i).z
		;simple_write_hdf5, nuv,gpstr + '/G_ABmag_NUV', fid
		;simple_write_hdf5, mu, gpstr + '/G_ABmag_u', fid
		;simple_write_hdf5, mg, gpstr + '/G_ABmag_g', fid
		;simple_write_hdf5, mr, gpstr + '/G_ABmag_r', fid
		;simple_write_hdf5, mi, gpstr + '/G_ABmag_i', fid
		;simple_write_hdf5, mz, gpstr + '/G_ABmag_z', fid
		

		FOR fi=0L, N_ELEMENTS(settings.flux_list)-1L DO BEGIN
			nullarr	= -1.d
			txtdum	= '/G_SB_' + STRTRIM(settings.flux_list(fi),2)

			IF N_ELEMENTS(bprop.abmag) GE 2L THEN void2	= EXECUTE('nullarr = bprop.sb(' + STRTRIM(i,2) + ').' + settings.flux_list(fi))
			
			simple_write_hdf5, nullarr, gpstr + txtdum, fid
		ENDFOR

		;nuv = -1.d & mu = -1.d & mg = -1.d & mr = -1.d & mi = -1.d & mz = -1.d
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN nuv= bprop.sb(i).nuv
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mu = bprop.sb(i).u
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mg = bprop.sb(i).g
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mr = bprop.sb(i).r
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mi = bprop.sb(i).i
		;IF N_ELEMENTS(bprop.ABmag) GE 2L THEN mz = bprop.sb(i).z
		;simple_write_hdf5, nuv, gpstr + '/G_SB_NUV', fid
		;simple_write_hdf5, mu,  gpstr + '/G_SB_u', fid
		;simple_write_hdf5, mg,  gpstr + '/G_SB_g', fid
		;simple_write_hdf5, mr,  gpstr + '/G_SB_r', fid
		;simple_write_hdf5, mi,  gpstr + '/G_SB_i', fid
		;simple_write_hdf5, mz,  gpstr + '/G_SB_z', fid

		cm = -1.d & cn = -1.d
		IF N_ELEMENTS(bprop.confrac_m) GE 2L THEN cm = REFORM(bprop.confrac_m(i).aper, nconf)
		IF N_ELEMENTS(bprop.confrac_n) GE 2L THEN cn = REFORM(bprop.confrac_n(i).aper, nconf)
		simple_write_hdf5, cm, gpstr + '/G_ConFrac_M', fid
		simple_write_hdf5, cn, gpstr + '/G_ConFrac_N', fid
		

		isclump = -1L
		IF N_ELEMENTS(bprop.isclump) GE 2L THEN isclump = bprop.isclump(i)
		simple_write_hdf5, isclump, idstr + '/isclump', fid
	
		;;----- Particle ID
		simple_write_hdf5, ptcl_id,	ppstr + '/P_ID',		fid

		;;----- Write Other properties
		;rate = -1.
		;IF N_ELEMENTS((*runstat.rv_ptmatch).rate) GE 2L THEN $
		;	rate = (*runstat.rv_ptmatch).rate(i)
		;simple_write_hdf5, rate, idstr + '/rate', fid

		simple_write_hdf5, (*runstat.rv_ptmatch).a_exp, idstr + '/Aexp', fid

		dom_list = -1L
		IF N_ELEMENTS((*runstat.rv_ptmatch).dom_list) GE 2L THEN $
			dom_list = REFORM((*runstat.rv_ptmatch).dom_list(i,*), nmpi)
		simple_write_hdf5, dom_list, idstr + '/Domain_List', fid

	ENDFOR
	H5F_CLOSE, fid
	;SPAWN, 'chmod 777 ' + STRTRIM(fname) + '/GAL_*.hdf5'
End

