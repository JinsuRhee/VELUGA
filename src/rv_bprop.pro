FUNCTION rv_BProp, settings, veluga, runstat, run=run

dir_data	= runstat.dir

;;-----
;; Check procedure set
;;-----
IF run EQ 0L THEN RETURN, PTR_NEW({ABMag:-1, SFR:-1, SFR_R:-1, SFR_T:-1, MAG_R:-1, isclump:-1., SFR_clumpcorr:-1., confrac:-1.},/no_copy)
IF run EQ 1L THEN BEGIN
	IF STRLEN(FILE_SEARCH(dir_data + 'rv_bprop.sav')) GE 5L THEN BEGIN
		RESTORE, dir_data + 'rv_bprop.sav'
		RETURN, PTR_NEW(output,/no_copy)
	ENDIF ELSE BEGIN
		run	= 2L
	ENDELSE
ENDIF
IF run EQ 2L THEN BEGIN
	PRINT, '        %%%%% (No previous works are found)'

	rawdata	= *runstat.rv_raw
	idlist	= *runstat.rv_id
	ptdata	= *runstat.rv_ptmatch
	;;-----
	;; Settings
	;;-----
	n_gal	= N_ELEMENTS(rawdata.id)
	n_part	= N_ELEMENTS(idlist.p_id)
	n_flux	= N_ELEMENTS(settings.flux_list)
	n_sfr	= N_ELEMENTS(settings.SFR_R)
	n_magap	= N_ELEMENTS(settings.MAG_R)

	;;-----
	;; Allocate Memory
	;;-----
	fl		= DBLARR(n_part, N_ELEMENTS(settings.flux_list)) - 1.0d8

	sfactor		= DBLARR(n_part)
	gyr		= DBLARR(n_part)

	abmag		= DBLARR(n_gal, n_flux, n_magap)
	SFR		= DBLARR(n_gal, n_sfr)

	;confrac		= DBLARR(n_gal, N_ELEMENTS(settings.CONF_r)) - 1.0d8
	PRINT, '        %%%%% BProp - MEMORY ALLOCATED'

	;;-----
	;; Conformal Time to SFactor and Gyr
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		dummy	= veluga->g_gyr(ptdata.p_age, runstat.snap)
		sfactor	= dummy.sfact
		gyr	= dummy.gyr
		PRINT, '        %%%%% BProp - CONFORMAL TIME CONVERTED'
	ENDIF

	;;-----
	;; SFR Calculation
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		mass 	= ptdata.p_mass

		FOR i=0L, n_gal-1L DO BEGIN
			FOR j=0L, n_sfr-1L DO BEGIN
				dtime 	= settings.SFR_T(i)
				dsize	= settings.SFR_R(i) * rawdata.r_halfmass(i)

				
				gdum	= [gyr(bind(i,0):bind(i,1), gyr(uind(i:0):uind(i,1)]
				mdum	= [mass(bind(i,0):bind(i,1), mass(uind(i:0):uind(i,1)]

				xdum	= [ptdata.p_pos(bind(i,0):bind(i,1),0):ptdata.p_pos(uind(i,0):uind(i,1),0)]
				ydum	= [ptdata.p_pos(bind(i,0):bind(i,1),1):ptdata.p_pos(uind(i,0):uind(i,1),1)]
				zdum	= [ptdata.p_pos(bind(i,0):bind(i,1),2):ptdata.p_pos(uind(i,0):uind(i,1),2)]

				sfr(i,j)	= veluga->g_sfr(xdum, ydum, zdum, gdum, mdum, rawdata.xc(i), rawdata.yc(i), rawdata.zc(i), aperture=dsize, timewindow=dtime)
			ENDFOR
		ENDFOR
	ENDIF

	;;FROM HERE
	;;-----
	;; CLUMP CORRECTION
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		cut	= WHERE(SFR(*,0) GT $
			rawdata.mass_tot / (settings.SFR_T(0)*1e9) * settings.clump_mfrac, ncut)
		SFR2	= SFR
		isclump	= LONARR(N_ELEMENTS(rawdata.id)) - 1L
		IF ncut GE 1L THEN BEGIN
			isclump(cut)	= 1L

			FOR i=0L, ncut-1L DO BEGIN
				ind	= cut(i)
				hostid	= rawdata.hostHaloID(ind)
				IF hostid LT 0L THEN CONTINUE
				cut2	= WHERE(rawdata.ID EQ hostid)
				SFR2(cut2,*)	+= SFR(ind,*)
				SFR2(ind,*)	= 0.
			ENDFOR
		ENDIF
		output	= CREATE_STRUCT('SFR', SFR, 'SFR_clumpcorr', SFR2, 'isclump', isclump)
		PRINT, '        %%%%% BProp - SFRs are calculated'
	ENDIF

	;;-----
	;; Magnitude
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		abmag	= get_mag(rawdata.xc, rawdata.yc, rawdata.zc, rawdata.r_halfmass, $
			idlist.b_ind, idlist.u_ind, ptdata.p_pos, ptdata.p_met, gyr, ptdata.p_mass, $
			MAG_R=settings.MAG_R, flux_list=settings.flux_list, $
			lib=settings.dir_lib, num_thread=settings.num_thread)

		output	= CREATE_STRUCT(output, 'ABMag', abmag)

		output	= CREATE_STRUCT(output, 'SFR_R', settings.SFR_R, 'SFR_T', settings.SFR_T, $
			'MAG_R', settings.MAG_R)
		PRINT, '        %%%%% BProp - Magnitudes are calculated'
	ENDIF
	;;-----
	;; Contamination Fraction
	;;-----
	confrac	= get_cfrac(settings, rawdata, n_snap, horg=settings.horg)

	output	= CREATE_STRUCT(output, 'CONFRAC_M', confrac.m)
	output	= CREATE_STRUCT(output, 'CONFRAC_N', confrac.n)
	PRINT, '        %%%%% BProp - Contamination fractions are calculated'

	SAVE, filename=dir_data + 'rv_gprop.sav', output
	RETURN, PTR_NEW(output,/no_copy)
ENDIF
END
