FUNCTION rv_ptmatch, settings, veluga, runstat, run=run

dir_data	= runstat.dir

;;-----
;; Check procedure set
;;-----
IF run EQ 0L THEN RETURN, PTR_NEW({p_pos:-1, p_vel:-1, p_age:-1, p_met:-1, p_mass:-1, $
	dom_list:-1, rate:-1, a_exp:-1},/no_copy)
IF run EQ 1L THEN BEGIN
	IF STRLEN(FILE_SEARCH(dir_data + '/rv_ptmatch.sav')) GE 5L THEN BEGIN
		RESTORE, dir_data + '/rv_ptmatch.sav'
		RETURN, PTR_NEW(output,/no_copy)
	ENDIF ELSE BEGIN
		run	= 2L
	ENDELSE
ENDIF
IF run EQ 2L THEN BEGIN
	;PRINT, '        %%%%% (No previous works are found)'

	output	= *runstat.rv_raw
	idlist	= *runstat.rv_id

	;;-----
	;; Read Simulation Info
	;;-----
	;infoname	= settings.dir_raw + 'output_' + STRING(runstat.snap,format='(I5.5)') + $
	;	'/info_' + STRING(runstat.snap,format='(I5.5)') + '.txt'
	siminfo 	= veluga->g_info(runstat.snap);, siminfo, file=infoname

	n_mpi	= N_ELEMENTS(siminfo.hindex(*,0))

        xc      = DBLARR(N_ELEMENTS(output.id))
        yc      = DBLARR(N_ELEMENTS(output.id))
        zc      = DBLARR(N_ELEMENTS(output.id))
        rr      = DBLARR(N_ELEMENTS(output.id))
        for i=0L, N_ELEMENTS(xc)-1L do xc(i) = output.xc(i)  * 3.086d21 / siminfo.unit_l
        for i=0L, N_ELEMENTS(yc)-1L do yc(i) = output.yc(i)  * 3.086d21 / siminfo.unit_l
        for i=0L, N_ELEMENTS(zc)-1L do zc(i) = output.zc(i)  * 3.086d21 / siminfo.unit_l
        for i=0L, N_ELEMENTS(rr)-1L do rr(i) = output.r_halfmass(i)* 3.086d21 / siminfo.unit_l
	rr	= rr*0.0 + MEDIAN(rr)

	rate	= FLTARR(N_ELEMENTS(output.id))
	n_bdn   = idlist.b_ind(-1,-1)+1L
	n_ubd   = idlist.u_ind(-1,-1) - idlist.b_ind(-1,-1)
	n_part  = n_bdn + n_ubd

	pos_pt  = DBLARR(n_part,3) - 1d8
	vel_pt  = pos_pt
	met_pt  = DBLARR(n_part) - 1d8
	age_pt  = met_pt
	mp_pt   = met_pt

	dom_list	= LONARR(N_ELEMENTS(xc),n_mpi) - 1L

	;;-----
	;; Matching Start
	;;-----
	n_nomatch	= N_ELEMENTS(WHERE(rate LT 0.8))
	match_cutval	= 0.9999
	N_itr = 0L & N_itrmax = 10L & dfact = 20.0

	dmp_mass	= 1./(DOUBLE(settings.neff)^3) * $
		(siminfo.omega_m - siminfo.omega_b) / siminfo.omega_m

	REPEAT BEGIN
		cut	= WHERE(rate LT match_cutval, ncut)

		PRINT, "%123123 ----------"
		PRINT, "%     ", STRTRIM(N_itr+1,2) + ' th iterations'
		PRINT, "%     FOR", ncut, " GALAXIES"
		PRINT, "%     USING ", dfact
		IF ncut EQ 0L THEN N_itr = N_itrmax


		IF N_itr LT N_itrmax THEN BEGIN
			;;----- Allocate Mem
			xc2 = xc(cut) & yc2 = yc(cut) & zc2 = zc(cut) & rr2 = rr(cut)

			rate2	= rate(cut)
			ind_b2 = idlist.b_ind(cut,*) & ind_u2 = idlist.u_ind(cut,*)

			id_pt2	= LON64ARR(1) + 0
			m0	= 0L
			FOR i2=0L, ncut - 1L DO BEGIN
				n0	= m0
				veluga->g_makearr, id_pt2, $
					[idlist.p_id(ind_b2(i2,0):ind_b2(i2,1))], m0, $
					unitsize=1000000L, type='L64'
				ind_b2(i2,0) = n0 & ind_b2(i2,1) = m0 - 1L
			ENDFOR

			FOR i2=0L, ncut - 1L DO BEGIN
				n0	= m0
				veluga->g_makearr, id_pt2, $
					[idlist.p_id(ind_u2(i2,0):ind_u2(i2,1))], m0, $
					unitsize=1000000L, type='L64'
				ind_u2(i2,0) = n0 & ind_u2(i2,1) = m0 - 1L
			ENDFOR
			id_pt2 = id_pt2(0L:m0-1L) & n_part2 = m0
			pos_pt2	= DBLARR(n_part2,3) - 1.0d8
			vel_pt2 = DBLARR(n_part2,3) - 1.0d8
			met_pt2 = DBLARR(n_part2) - 1.0d8
			age_pt2 = met_pt2 & mp_pt2 = met_pt2

			;;----- Domain
			ftr_name	= settings.dir_lib + 'src/fortran/find_domain.so'
				dom_list2	= LONARR(N_ELEMENTS(xc2),n_mpi) - 1L
				larr	= LONARR(20) & darr = DBLARR(20)
				larr(0) = N_ELEMENTS(xc2)
				larr(1) = n_mpi
				larr(2)	= 1L;settings.num_thread

				;IF dfact GE 100. THEN $
				;	darr(0) = dfact / (max(rr2) * siminfo.unit_l / 3.086d21)
				;	;; Search in Physical Radius

				darr(0) = dfact		;; Radius factor

			void	= CALL_EXTERNAL(ftr_name, 'find_domain', $
				xc2, yc2, zc2, rr2, siminfo.hindex, siminfo.levmax, $
				dom_list2, larr, darr)

			;;----- Matching
			ftr_name	= settings.dir_lib + 'src/fortran/rv_match.so'
			lset = LONARR(20) & dset = DBLARR(20)
				lset(0) = N_ELEMENTS(id_pt2)	;; # of particles in VR Data
				lset(1) = N_ELEMENTS(xc2)	;; # of Gals
				lset(2) = settings.num_thread
				lset(3) = 10L			;; # of Doamins in a set
				lset(4) = n_mpi
				lset(5) = runstat.snap
				lset(10) = STRLEN(settings.dir_raw)
				IF settings.famtype EQ 'old' THEN lset(18) = 100 ;; For YZiCS ver
				IF settings.idtype EQ 'long' THEN lset(19) = 100 ;; For logn int ID

				dset(0) = dmp_mass
			IF settings.horg EQ 'g' THEN dset(1) = 1.
			IF settings.horg EQ 'h' THEN dset(1) = -1.

			void = CALL_EXTERNAL(ftr_name, 'rv_match', $
				lset, dset, settings.dir_raw, $
				id_pt2, ind_b2, ind_u2, pos_pt2, vel_pt2, $
				met_pt2, age_pt2, mp_pt2, rate2, $
				dom_list2)

			;;----- Extract
			FOR i2=0L, N_ELEMENTS(cut) - 1L DO BEGIN
				IF N_itr LT N_itrmax - 1L THEN $
					IF rate2(i2) LT match_cutval THEN CONTINUE

				rate(cut(i2)) = rate2(i2)

				pos_pt(idlist.b_ind(cut(i2),0):idlist.b_ind(cut(i2),1),*) = $
					pos_pt2(ind_b2(i2,0):ind_b2(i2,1),*)
				vel_pt(idlist.b_ind(cut(i2),0):idlist.b_ind(cut(i2),1),*) = $
					vel_pt2(ind_b2(i2,0):ind_b2(i2,1),*)
				met_pt(idlist.b_ind(cut(i2),0):idlist.b_ind(cut(i2),1)) = $
					met_pt2(ind_b2(i2,0):ind_b2(i2,1))
				age_pt(idlist.b_ind(cut(i2),0):idlist.b_ind(cut(i2),1)) = $
					age_pt2(ind_b2(i2,0):ind_b2(i2,1))
				mp_pt(idlist.b_ind(cut(i2),0):idlist.b_ind(cut(i2),1)) = $
					mp_pt2(ind_b2(i2,0):ind_b2(i2,1))

				pos_pt(idlist.u_ind(cut(i2),0):idlist.u_ind(cut(i2),1),*) = $
					pos_pt2(ind_u2(i2,0):ind_u2(i2,1),*)
				vel_pt(idlist.u_ind(cut(i2),0):idlist.u_ind(cut(i2),1),*) = $
					vel_pt2(ind_u2(i2,0):ind_u2(i2,1),*)
				met_pt(idlist.u_ind(cut(i2),0):idlist.u_ind(cut(i2),1)) = $
					met_pt2(ind_u2(i2,0):ind_u2(i2,1))
				age_pt(idlist.u_ind(cut(i2),0):idlist.u_ind(cut(i2),1)) = $
					age_pt2(ind_u2(i2,0):ind_u2(i2,1))
				mp_pt(idlist.u_ind(cut(i2),0):idlist.u_ind(cut(i2),1)) = $
					mp_pt2(ind_u2(i2,0):ind_u2(i2,1))

				dom_list(cut(i2),*) = dom_list2(i2,*)
			ENDFOR

			N_itr ++
			IF dfact LT 20.0 THEN BEGIN
				dfact = dfact + 40.0
			ENDIF ELSE BEGIN
				dfact = dfact + 50.0
			ENDELSE
		ENDIF

		void	= WHERE(rate LT match_cutval, leftgal)
		PRINT, "%     Done. Left Galaxies are ", leftgal
		PRINT, "%123123 ----------"

		IF leftgal EQ 0L THEN N_itr = N_itrmax
	ENDREP UNTIL N_itr GE N_itrmax

	pos_pt	= pos_pt * siminfo.unit_l / 3.086d21
	vel_pt	= vel_pt * siminfo.kms
	mp_pt	= mp_pt * siminfo.unit_m / 1.98892e33

	output	= CREATE_STRUCT('p_pos', pos_pt, $
		'p_vel', vel_pt, 'p_age', age_pt, $
		'p_met', met_pt, 'p_mass', mp_pt, $
		'dom_list', dom_list, 'rate', rate, $
		'a_exp', siminfo.aexp)

	IF settings.pp_saveprocess EQ 1L THEN SAVE, filename=dir_data + '/rv_ptmatch.sav', output
	RETURN, PTR_NEW(output,/no_copy)
ENDIF
END
