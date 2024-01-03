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
	;PRINT, '        %%%%% (No previous works are found)'

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
	n_aper 	= N_ELEMENTS(settings.CONF_R)
	;;-----
	;; Allocate Memory
	;;-----
	fl		= DBLARR(n_part, N_ELEMENTS(settings.flux_list)) - 1.0d8

	sfactor	= DBLARR(n_part)
	gyr		= DBLARR(n_part)

	abmag0	= '{'
	FOR i=0L, n_flux-1L DO BEGIN
		abmag0 = abmag0 + settings.flux_list(i) + ':DBLARR(' + STRING(n_magap) + ')'
		IF i NE n_flux-1L THEN abmag0 = abmag0 + ', '
	ENDFOR
	abmag0	= abmag0 + '}'
	void 	= EXECUTE('abmag0 = ' + abmag0)
	void 	= EXECUTE('abmag = REPLICATE(abmag0, ' +  STRING(n_gal) + ')')

	SFR		= DBLARR(n_gal, n_sfr)

	
	aperdum = DBLARR(n_gal, n_aper)
	;confrac	= REPLICATE(n_gal, {n:0.d, m:0.d}) - 1.0d8
	
	;PRINT, '        %%%%% BProp - MEMORY ALLOCATED'

	;;-----
	;; Conformal Time to SFactor and Gyr
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		dummy	= veluga->g_gyr(runstat.snap, ptdata.p_age)
		sfactor	= dummy.sfact
		gyr		= dummy.gyr
		metal 	= ptdata.p_met

		veluga->ppout2, ' - Conformal time converted'
		;PRINT, '        %%%%% BProp - CONFORMAL TIME CONVERTED'
	ENDIF

	;;-----
	;; SFR Calculation
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		mass 	= ptdata.p_mass

		FOR i=0L, n_gal-1L DO BEGIN
			FOR j=0L, n_sfr-1L DO BEGIN
				dtime 	= settings.SFR_T(j)
				dsize	= settings.SFR_R(j) * rawdata.r_halfmass(i)

				
				gdum	= [gyr(idlist.b_ind(i,0):idlist.b_ind(i,1)), gyr(idlist.u_ind(i,0):idlist.u_ind(i,1))]
				mdum	= [mass(idlist.b_ind(i,0):idlist.b_ind(i,1)), mass(idlist.u_ind(i,0):idlist.u_ind(i,1))]

				xdum	= [ptdata.p_pos(idlist.b_ind(i,0):idlist.b_ind(i,1),0), ptdata.p_pos(idlist.u_ind(i,0):idlist.u_ind(i,1),0)]
				ydum	= [ptdata.p_pos(idlist.b_ind(i,0):idlist.b_ind(i,1),1), ptdata.p_pos(idlist.u_ind(i,0):idlist.u_ind(i,1),1)]
				zdum	= [ptdata.p_pos(idlist.b_ind(i,0):idlist.b_ind(i,1),2), ptdata.p_pos(idlist.u_ind(i,0):idlist.u_ind(i,1),2)]

				sfr(i,j)	= veluga->g_sfr(xdum, ydum, zdum, gdum, mdum, rawdata.xc(i), rawdata.yc(i), rawdata.zc(i), aperture=dsize, timewindow=dtime)

			ENDFOR
		ENDFOR

		veluga->ppout2, ' - SFR computed'
	ENDIF


	;;-----
	;; CLUMP CORRECTION
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		cut	= WHERE(SFR(*,0) GT $
			rawdata.mass_tot / (settings.SFR_T(0)*1e9) * settings.pp_clump_mfrac, ncut)
		SFR2	= SFR
		isclump	= LONARR(N_ELEMENTS(rawdata.id)) - 1L
		IF ncut GE 1L THEN BEGIN
			isclump(cut)	= 1L

		;	FOR i=0L, ncut-1L DO BEGIN
		;		ind	= cut(i)
		;		hostid	= rawdata.hostHaloID(ind)
		;		IF hostid LT 0L THEN CONTINUE
		;		cut2	= WHERE(rawdata.ID EQ hostid)
		;		SFR2(cut2,*)	+= SFR(ind,*)
		;		SFR2(ind,*)	= 0.
		;	ENDFOR
		ENDIF
		output	= CREATE_STRUCT('SFR', SFR, 'SFR_clumpcorr', SFR2, 'isclump', isclump)
		veluga->ppout2, ' - Clump classification computed'
	ENDIF

	;;-----
	;; Magnitude
	;;-----
	IF settings.horg EQ 'g' THEN BEGIN
		FOR gi=0L, n_gal-1L DO BEGIN
			gdum	= [gyr(idlist.b_ind(gi,0):idlist.b_ind(gi,1)), gyr(idlist.u_ind(gi,0):idlist.u_ind(gi,1))]
			mdum	= [mass(idlist.b_ind(gi,0):idlist.b_ind(gi,1)), mass(idlist.u_ind(gi,0):idlist.u_ind(gi,1))]
			metdum 	= [metal(idlist.b_ind(gi,0):idlist.b_ind(gi,1)), metal(idlist.u_ind(gi,0):idlist.u_ind(gi,1))]

			xdum	= [ptdata.p_pos(idlist.b_ind(gi,0):idlist.b_ind(gi,1),0), ptdata.p_pos(idlist.u_ind(gi,0):idlist.u_ind(gi,1),0)]
			ydum	= [ptdata.p_pos(idlist.b_ind(gi,0):idlist.b_ind(gi,1),1), ptdata.p_pos(idlist.u_ind(gi,0):idlist.u_ind(gi,1),1)]
			zdum	= [ptdata.p_pos(idlist.b_ind(gi,0):idlist.b_ind(gi,1),2), ptdata.p_pos(idlist.u_ind(gi,0):idlist.u_ind(gi,1),2)]

			d3d 	= veluga->g_d3d(xdum, ydum, zdum, [rawdata.xc(gi), rawdata.yc(gi), rawdata.zc(gi)])

			ldum 	= veluga->g_luminosity(mdum, gdum, metdum, settings.flux_list)

			FOR ai=0L, n_magap-1L DO BEGIN
				dsize 	= settings.MAG_R(ai) * rawdata.r_halfmass(gi)

				IF dsize LT 0. THEN BEGIN
					ldum2 	= ldum
				ENDIF ELSE BEGIN
					cut 	= WHERE(d3d LT dsize, ncut)
					ldum2 	= ldum(cut)
				ENDELSE

				
				FOR fi=0L, n_flux-1L DO BEGIN
					CASE settings.flux_list(fi) OF
						'u' 	: abmag(gi).u(ai) 	= veluga->g_mag(ldum.u, ['u'])
						'g' 	: abmag(gi).g(ai) 	= veluga->g_mag(ldum.g, ['g'])
						'r' 	: abmag(gi).r(ai) 	= veluga->g_mag(ldum.r, ['r'])
						'i' 	: abmag(gi).i(ai) 	= veluga->g_mag(ldum.i, ['i'])
						'z' 	: abmag(gi).z(ai) 	= veluga->g_mag(ldum.z, ['z'])
						'nuv' 	: abmag(gi).nuv(ai) = veluga->g_mag(ldum.nuv, ['nuv'])
						'NUV' 	: abmag(gi).nuv(ai) = veluga->g_mag(ldum.nuv, ['nuv'])
					ENDCASE
				ENDFOR
			ENDFOR
		ENDFOR

		output 	= CREATE_STRUCT(output, 'ABMag', abmag)
		veluga->ppout2, ' - Magnitude computed'
	ENDIF

	;;-----
	;; Contamination Fraction
	;;-----

	IF settings.horg EQ 'g' THEN BEGIN
		sizedum 	= rawdata.r_halfmass
	ENDIF ELSE BEGIN
		sizedum 	= rawdata.rvir
	ENDELSE

	FOR i=0L, n_aper-1L DO aperdum(*,i) = sizedum * settings.CONF_r(i)
	
	cdum 	= veluga->g_cfrac(runstat.snap, rawdata.xc, rawdata.yc, rawdata.zc, aperdum)
	;cdum 	= veluga->g_cfrac(runstat.snap, rawdata.xc(0:10), rawdata.yc(0:10), rawdata.zc(0:10), sizedum(0:10), aperdum(0:10,*))

	cdum_n	= REPLICATE({aper:DBLARR(n_aper)}, n_gal)
	cdum_m	= REPLICATE({aper:DBLARR(n_aper)}, n_gal)
	FOR i=0L, n_aper-1L DO cdum_n.aper(i) = cdum.n(*,i)
	FOR i=0L, n_aper-1L DO cdum_m.aper(i) = cdum.m(*,i)
	output	= CREATE_STRUCT(output, 'CONFRAC_M', cdum_n)
	output	= CREATE_STRUCT(output, 'CONFRAC_N', cdum_m)
	veluga->ppout2, ' - Contamination Fraction computed'
	
	STOP

	SAVE, filename=dir_data + 'rv_gprop.sav', output
	RETURN, PTR_NEW(output,/no_copy)
ENDIF
END
