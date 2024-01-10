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

