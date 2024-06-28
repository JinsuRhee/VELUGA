FUNCTION veluga::init, fname, num_thread=num_thread

	IF ~KEYWORD_SET(num_thread) THEN num_thread = 1L

	self.header	= PTR_NEW(/allocate)
	self.num_thread	= num_thread


	self.rdheader, fname
	ptr 	= REPLICATE({tree:PTR_NEW(), key:PTR_NEW(), stat:-2L},2)
	self.tree 	= PTR_NEW(ptr)
	
	settings 	= self.getheader()

	RETURN, 1
END

FUNCTION veluga::allocate, nn, type=type

	settings 	= self->getheader()
	CASE type OF
		'part'		: BEGIN
			;RETURN, REPLICATE({xx:0.d, yy:0.d, zz:0.d, vx:0.d, vy:0.d, vz:0.d, mp:0.d, ap:0.d, zp:0.d, gyr:0.d, redsh:0.d, sfact:0.d, id:0L, family:0L, domain:0L, KE:0.d, UE:0.d, PE:0.d}, nn)
			da 	= DBLARR(nn)
			la 	= LONARR(nn)
			pa 	= PTRARR(nn)
			RETURN, {N:nn, xx:da, yy:da, zz:da, vx:da, vy:da, vz:da, mp:da, ap:da, zp:da, gyr:da, redsh:da, sfact:da, id:la, family:la, domain:la, KE:da, UE:da, PE:da, dum1:pa, dum2:pa, dum3:pa, dum4:pa, dum5:pa, tag:'part'}
			END
		'cell'		: BEGIN
			da 	= DBLARR(nn)
			la 	= LONARR(nn)
			pa 	= PTRARR(nn)
			info 	= self->g_info(1L)

			IF N_ELEMENTS(settings.hydro_variables) LE 7L THEN BEGIN ;; specify by the # of elements
				RETURN, {N:nn, xx:da, yy:da, zz:da, vx:da, vy:da, vz:da, level:la, dx:da, den:da, temp:da, zp:da, mp:da, KE:da, UE:da, PE:da, P_thermal:da, levelind:LONARR(info.levmax+1L,3), dum1:pa, dum2:pa, dum3:pa, dum4:pa, dum5:pa, tag:'cell'}
			ENDIF ELSE BEGIN
					additional_hvar_tag 	= settings.hydro_variables(6L:*)
					tmp 	= {N:nn, xx:da, yy:da, zz:da, vx:da, vy:da, vz:da, level:la, dx:da, den:da, temp:da, zp:da, mp:da, KE:da, UE:da, PE:da, P_thermal:da, levelind:LONARR(info.levmax+1L,3), dum1:pa, dum2:pa, dum3:pa, dum4:pa, dum5:pa, tag:'cell'}
					FOR i=0L, N_ELEMENTS(additional_hvar_tag)-1L DO $
						tmp 	= CREATE_STRUCT(tmp, additional_hvar_tag(i), da)

					RETURN, tmp;REPLICATE(tmp, nn)
				ENDELSE
		    END
		;'cell_eff' 	: BEGIN

		ELSE: STOP
	ENDCASE
END
;;-----
;; HEADER
;;-----
PRO veluga::rdheader, fname

	FINDPRO, 'veluga__define', dirlist=dirlist, /noprint
	
	settings 	= {header:'^^', dir_lib:dirlist(0), num_thread:1L, sun_met:0.d, $
		erg_msg_i:'%123123----- VELUGA -----', erg_msg_f:'%123123----------', erg_msg_0:' ', $
		pp_msg_i:'%123123----- VELUGA (Post processing) -----', pp_msg_f:'%123123----------', pp_msg_0:' '}

	
	OPENR, 10, fname
	FOR i=0L, FILE_LINES(fname)-1L DO BEGIN
		v1 	= STRING(' ')
		READF, 10, v1
		v1 	= STRTRIM(v1, 2)

		IF STRLEN(v1) EQ 0L THEN CONTINUE
		in      = STRPOS(v1, '#')
		IF MAX(in) GE 0L THEN CONTINUE

		void	= EXECUTE(v1)
		tag_name= (STRSPLIT(v1, '=', /extract))[0]
		v2      = 'settings = CREATE_STRUCT(settings, "' + STRTRIM(tag_name,2) + '", ' + $
			STRTRIM(tag_name,2) + ')'
		void	= EXECUTE(v2)
	ENDFOR
	CLOSE, 10
	settings.num_thread	= self.num_thread
	settings.sun_met	= 0.02d

	(*self.header)	= settings
END


FUNCTION veluga::getheader
	RETURN, *self.header
END



PRO veluga::errorout, str
	settings 	= self->getheader()
	PRINT, settings.erg_msg_i
	PRINT, settings.erg_msg_0
	PRINT, str
	PRINT, settings.erg_msg_0
	PRINT, settings.erg_msg_f
	RETURN
END

PRO veluga::ppout, str
	settings 	= self->getheader()
	PRINT, settings.pp_msg_i
	PRINT, settings.pp_msg_0
	PRINT, str
	PRINT, settings.pp_msg_0
	PRINT, settings.pp_msg_f
	RETURN
END

PRO veluga::ppout2, str
	settings 	= self->getheader()
	PRINT, '		', str
	RETURN
END

;;-----
;; READ ROUTINE
;;-----
FUNCTION veluga::r_gal_getdata, fid, str
	did 	= H5D_OPEN(fid, str)
	dumarr 	= H5D_READ(did)
	H5D_CLOSE, did
	RETURN, dumarr
END
FUNCTION veluga::r_gal, snap0, id0, horg=horg, Gprop=Gprop
	;;-----
	;;
	;; IF Grpop is argued, only selected field are read
	;;-----

	;;-----
	;; READ Galaxies
	;;-----
	IF ~KEYWORD_SET(horg) THEN horg = 'g'
	settings	= self.getheader()
	dir 		= settings.dir_catalog

	IF horg EQ 'g' THEN $
		fname = dir + 'Galaxy/VR_Galaxy/snap_' + STRING(snap0,format='(I4.4)') + '.hdf5'
	IF horg EQ 'h' THEN $
		fname = dir + 'Halo/VR_Halo/snap_' + STRING(snap0,format='(I4.4)') + '.hdf5'

	;;-----
	;; OPEN HDF5
	;;-----
	fid	= H5F_OPEN(fname)

	;;-----
	;; READ OVERALL INFO FIRST
	;;-----
	flux_list 	= self->r_gal_getdata(fid, 'Flux_List')
	sfr_r 		= self->r_gal_getdata(fid, 'SFR_R')
	sfr_t 		= self->r_gal_getdata(fid, 'SFR_T')
	mag_r 		= self->r_gal_getdata(fid, 'MAG_R')
	conf_r 		= self->r_gal_getdata(fid, 'CONF_R')

	mass_tot 	= self->r_gal_getdata(fid, 'Mass_tot')
	r_halfmass	= self->r_gal_getdata(fid, 'R_HalfMass')
	mvir 		= self->r_gal_getdata(fid, 'Mvir')
	rvir 		= self->r_gal_getdata(fid, 'Rvir')
	ID 			= self->r_gal_getdata(fid, 'ID')
	conf_m 		= self->r_gal_getdata(fid, 'CONF_M')
	conf_n 		= self->r_gal_getdata(fid, 'CONF_N')

	;;-----
	;; LOAD COLUMN
	;;-----
	IF ~KEYWORD_SET(gprop) THEN BEGIN 
		gprop		= settings.column_list;[settings.column_list, settings.gal_prop]
	
		FOR i=0L, N_ELEMENTS(settings.gal_prop)-1L DO BEGIN
	
			IF settings.gal_prop(i) EQ 'sfr' THEN settings.gal_prop(i) = 'SFR'
			IF settings.gal_prop(i) EQ 'SFR' THEN gprop	= [gprop, settings.gal_prop(i)]
	
			IF settings.gal_prop(i) EQ 'abmag' THEN settings.gal_prop(i) = 'ABmag'
			IF settings.gal_prop(i) EQ 'ABmag' THEN $
				FOR fi=0L, N_ELEMENTS(flux_list)-1L DO $
					gprop 	= [gprop, settings.gal_prop(i) + '_' + STRTRIM(flux_list(fi))]
			
	
			IF settings.gal_prop(i) EQ 'sb' THEN settings.gal_prop(i) = 'SB'
			IF settings.gal_prop(i) EQ 'SB' THEN $
				FOR fi=0L, N_ELEMENTS(flux_list)-1L DO $
					gprop 	= [gprop, settings.gal_prop(i) + '_' + STRTRIM(flux_list(fi))]
			
			IF settings.gal_prop(i) EQ 'conf' THEN settings.gal_prop(i) = 'CONF'
			IF settings.gal_prop(i) EQ 'confrac' THEN settings.gal_prop(i) = 'CONF'
			IF settings.gal_prop(i) EQ 'CONF' THEN $
				gprop 	= [gprop, 'ConFrac_M', 'ConFrac_N']
		ENDFOR
	ENDIF


	;;-----
	;; Initial cut (Not implemented Yet)
	;;-----
	IF id0 LT 0L THEN BEGIN 
	;		IF KEYWORD_SET(masscut) THEN BEGIN
	;			IF horg EQ 'g' THEN $
	;				mcut	= WHERE(mass_tot GE masscut(0) AND mass_tot LT masscut(1), nn)
	;			IF horg EQ 'h' THEN $
	;				mcut	= WHERE(mvir GE masscut(0) AND mvir LT masscut(1), nn)
	;
	;			IF nn EQ 0L THEN BEGIN
	;				PRINT, '%123123-----'
	;				PRINT, '	NO GALAXIES (HALOS) WITHIN THE MASS LIMIT'
	;				PRINT, '	(CONVERT TO READING ALL GALS'
	;
	;				mcut	= WHERE(mvir GE 0. AND mvir LT 1e20, nn)
	;			ENDIF
	;			n_gal	= nn
	;		ENDIF ELSE BEGIN
			mcut	= WHERE(mvir GE 0. AND mvir LT 1e20, nn)
			n_gal	= nn
	;		ENDELSE
	ENDIF ELSE BEGIN
		mcut	= WHERE(ID EQ id0, nn)
		n_gal	= nn
	ENDELSE

	;;-----
	;; Allocate Memory
	;;-----
	tmpstr	= 'GP = {snapnum:0L, redsh:0.d, '
	FOR j=0L, N_ELEMENTS(Gprop) - 1L DO BEGIN
		did	= H5D_OPEN(fid, 'ID_000001/G_Prop/G_' + $
			STRTRIM(Gprop(j),2))

		tmp	= 't' + Gprop(j) + ' = H5D_READ(did)'
		void	= EXECUTE(tmp)

		tmp2	= 'n_t = N_ELEMENTS(t' + Gprop(j) + ')'
		void	= EXECUTE(tmp2)

		IF n_t EQ 1L THEN void = EXECUTE('t' + Gprop(j) + ' = t' + Gprop(j) + '(0)')
		H5D_CLOSE, did
		tmpstr	+= Gprop(j) + ':t'+ Gprop(j) + ', '
	ENDFOR

	dlist 	= self->r_gal_getdata(fid, 'ID_000001/Domain_List')
	tmpstr	+= 'isclump:-1L, Aexp:1.0d, Domain_List:dlist,  '
	n_mpi	= N_ELEMENTS(dlist)

	tmpstr	+= 'flux_List:flux_list, CONF_R:CONF_R, MAG_R:MAG_R, SFR_R:SFR_R, SFR_T:SFR_T}'

	void	= EXECUTE(tmpstr)
	GP	= REPLICATE(GP, n_gal)
	
	;;-----
	;; Read Galaxies
	;;-----
	FOR i=0L, n_gal - 1L DO BEGIN
		ind	= mcut(i)
		idstr	= 'ID_' + STRING(ID(ind), format='(I6.6)')
		FOR j=0L, N_ELEMENTS(Gprop) - 1L DO BEGIN
			did	= H5D_OPEN(fid, idstr + '/G_Prop/G_' + $
				STRTRIM(Gprop(j),2))
			tmp	= 'GP(' + STRTRIM(i,2) + ').' + STRTRIM(Gprop(j),2) + $
				'=H5D_READ(did)'
			void	= EXECUTE(tmp)
			H5D_CLOSE, did
		ENDFOR

		
		GP(i).isclump	= self->r_gal_getdata(fid, idstr + '/isclump')
		GP(i).Domain_list 	= self->r_gal_getdata(fid, idstr + '/Domain_List')
		GP(i).aexp 		= self->r_gal_getdata(fid, idstr + '/Aexp')
		GP(i).snapnum 	= snap0
		GP(i).redsh 	= 1./GP(i).aexp - 1.d		

		;did	= H5D_OPEN(fid, idstr + '/rate')
		;GP(i).rate	= H5D_READ(did)
		;H5D_CLOSE, did
	ENDFOR

	H5F_CLOSE, fid
	RETURN, GP
END

FUNCTION veluga::r_part, snap0, id0, horg=horg, simout=simout

	;;-----
	;; READ MEMBER Part
	;;	/simout 	- output as the simulation raw unit
	;;-----
	IF ~KEYWORD_SET(horg) THEN horg='g'
	settings 	= self->getheader()

	;;----- READ PTCL IDs & domain list
	pid 	= self->r_pid(snap0, id0, horg=horg)
	dom_list= self->r_domain(snap0, id0, horg=horg)

	

	;;----- READ PTCL
	n_ptcl	= N_ELEMENTS(pid)
	info 	= self->g_info(snap0)
	dmp_mass 	= 1./(settings.neff*1.d)^3 * (info.omega_M - info.omega_b) / info.omega_m
	pinfo	= DBLARR(n_ptcl,9) -  1.0d8

	ftr_name 	= settings.dir_lib + 'src/fortran/get_ptcl.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0) = n_ptcl
		larr(1) = N_ELEMENTS(dom_list)
		larr(2) = snap0
		larr(3)	= self.num_thread
		larr(10)= STRLEN(settings.dir_raw)
		IF horg EQ 'g' THEN larr(11) = 10L
		IF horg EQ 'h' THEN larr(11) = -10L
		larr(18)= 0L
		larr(19)= 0L

		IF settings.famtype EQ 'new' THEN larr(18) = 100L
		IF settings.idtype EQ 'long64' THEN larr(19)= 100L
		IF horg EQ 'h' THEN darr(11) = dmp_mass

		void	= CALL_EXTERNAL(ftr_name, 'get_ptcl', $
			larr, darr, settings.dir_raw, pid, pinfo, dom_list)

	

	;; POST PROCESSING
	cut	= WHERE(pinfo(*,0) gt -1.0d7, ncut)
	IF MAX(cut) LT 0 THEN BEGIN
		self->errorout, '		r_part: NO MATCHED PTCLs'
		RETURN, -1.
	ENDIF

	;output	= {rate: (ncut * 1.d ) / n_ptcl}

	output 	= self->allocate(ncut, type='part')
	;REPLICATE({xx:0.d, yy:0.d, zz:0.d, vx:0.d, vy:0.d, vz:0.d, mp:0.d, ap:0.d, zp:0.d, id:0L}, ncut)

	output.xx 	= pinfo(cut,0)
	output.yy 	= pinfo(cut,1)
	output.zz 	= pinfo(cut,2)

	output.vx 	= pinfo(cut,3)
	output.vy 	= pinfo(cut,4)
	output.vz 	= pinfo(cut,5)

	output.mp 	= pinfo(cut,6)
	output.ap 	= pinfo(cut,7)
	output.zp 	= pinfo(cut,8)

	output.id 	= pid(cut)

	IF ~KEYWORD_SET(simout) THEN BEGIN
		output.xx 	*= (info.unit_l/info.cgs.kpc)
		output.yy 	*= (info.unit_l/info.cgs.kpc)
		output.zz 	*= (info.unit_l/info.cgs.kpc)

		output.vx 	*= (info.kms)
		output.vy 	*= (info.kms)
		output.vz 	*= (info.kms)

		output.mp 	*= (info.unit_m / info.cgs.m_sun)
	ENDIF

	RETURN, output
END

FUNCTION veluga::r_pid, snap0, id0, horg=horg
	;;-----
	;; READ PTCL ID
	;;-----

	IF ~KEYWORD_SET(horg) THEN horg='g'
	settings 	= self->getheader()

	dir 	= settings.dir_catalog

	IF horg EQ 'h' THEN $
		dir	= dir + 'Halo/VR_Halo/'

	IF horg EQ 'g' THEN $
		dir	= dir + 'Galaxy/VR_Galaxy/'
	
	fname	= dir + 'snap_' + STRING(snap0,format='(I4.4)') + '.hdf5'

	fid = H5F_OPEN(fname) & did = H5D_OPEN(fid, 'ID_' + STRING(id0,format='(I6.6)') + '/P_Prop/P_ID')
	pid = H5D_READ(did) & H5D_CLOSE, did & H5F_CLOSE, fid

	RETURN, pid
END

FUNCTION veluga::r_domain, snap0, id0, horg=horg

	IF ~KEYWORD_SET(horg) THEN horg='g'
	settings 	= self->getheader()

	dir 	= settings.dir_catalog

	IF horg EQ 'h' THEN $
		dir	= dir + 'Halo/VR_Halo/'

	IF horg EQ 'g' THEN $
		dir	= dir + 'Galaxy/VR_Galaxy/'
	
	fname	= dir + 'snap_' + STRING(snap0,format='(I4.4)') + '.hdf5'

	fid = H5F_OPEN(fname) & did = H5D_OPEN(fid, 'ID_' + STRING(id0,format='(I6.6)') + '/Domain_List')
	dom_list = H5D_READ(did) & H5D_CLOSE, did & H5F_CLOSE, fid

	cut 	= WHERE(dom_list GE 1L) + 1L
	RETURN, cut
END

FUNCTION veluga::r_cell, snap0, id0, rsize, horg=horg
	;;-----
	;; READ AMR CELL around a galaxy
	;;-----

	gal 	= self->r_gal(snap0, id0, horg=horg)
	RETURN, self->g_cell(gal.xc, gal.yc, gal.zc, rsize, snap0)
END

FUNCTION veluga::r_tree_load, horg=horg

	IF horg EQ 'g' THEN tind = 0L
	IF horg EQ 'h' THEN tind = 1L
	;IF horg EQ 'g' THEN 

	IF (*self.tree)(tind).stat EQ -2L THEN BEGIN
		settings 	= self->getheader()

		IF horg EQ 'g' THEN fname = settings.dir_catalog + 'Galaxy/tree/ctree.sav'
		IF horg EQ 'h' THEN fname = settings.dir_catalog + 'Halo/tree/ctree.sav'

		isfile 	= FILE_SEARCH(fname)
		
		IF STRLEN(isfile) GE 5L THEN BEGIN
			RESTORE, fname

			(*self.tree)(tind).key 	= PTR_NEW(tree_key)
			(*self.tree)(tind).tree = PTR_NEW(complete_tree)
			(*self.tree)(tind).stat = 1L
			
		ENDIF ELSE BEGIN
			(*self.tree)(tind).stat = -1L
		ENDELSE

	ENDIF

	RETURN, (*self.tree)(tind)
	
END
FUNCTION veluga::r_tree, snap0, id0, horg=horg

	;;-----
	;; READ TREE
	;;-----
	tree 	= self->r_tree_load(horg=horg)

	IF tree.stat EQ -1L THEN BEGIN
		self->errorout, 'No tree data exists for ' + horg
		RETURN, 1L
	ENDIF
	
	key 	= (*tree.key)(0)
	keyval 	= snap0 + id0*key
	ind 	= (*tree.key)(keyval)
	IF ind EQ -1L THEN RETURN, 1L
	RETURN, *(*tree.tree)(ind)
END

FUNCTION veluga::r_evol, snap0, id0, horg=horg

	;;-----
	;; READ TREE
	;;-----

	tree 	= self->r_tree(snap0, id0, horg=horg)
	
	IF TYPENAME(tree) EQ 'LONG' THEN BEGIN
		self->errorout, 'No tree data exists for ' + horg
		RETURN, 1L
	ENDIF

	;;-----
	;; First Galaxy
	;;-----
	g0 	= self->r_gal(snap0, id0, horg=horg)

	n_gal 	= N_ELEMENTS(tree.id)
	g 	= REPLICATE(g0, n_gal)

	FOR i=0L, n_gal-1L DO BEGIN
		g(i) 	= self->r_gal(tree.snap(i), tree.id(i), horg=horg)
	ENDFOR

	RETURN, g
END

;;-----
;; SIMPLE GET FTNS
;;-----
FUNCTION veluga::g_wmean, xx, ww
	RETURN, TOTAL(xx * ww) / TOTAL(ww)
END

FUNCTION veluga::g_wstddev, xx, ww
	m 	= self->g_wmean(xx, ww)
	RETURN, SQRT( TOTAL( ww * (xx - m)^2 ) / TOTAL(ww) )
END

FUNCTION veluga::g_unique, array
	dum	= array & dum = dum(SORT(dum)) & dum = dum(UNIQ(dum)) & RETURN, dum
END

FUNCTION veluga::g_d3d, x, y, z, cen
	RETURN, SQRT((x-cen(0))^2 + (y-cen(1))^2 + (z-cen(2))^2)
END

FUNCTION veluga::g_cdist, zarr

	settings 	= self->getheader()
	info 	= self->g_info(1L)

	tbl 	= self->t_comv_load(info.omega_M, info.omega_L, info.H0)

	cdist 	= INTERPOL(tbl.comv, tbl.redsh, zarr)

	RETURN, cdist
END

PRO veluga::g_makearr, array, input, index, unitsize=unitsize, type=type
	IF ~KEYWORD_SET(unitsize) THEN unitsize = 10000L
	IF ~KEYWORD_SET(type) THEN type = 'F'

	IF type NE 'L' AND type NE 'D' AND type NE 'F' AND type NE 'L64' THEN BEGIN
		PRINT, '***** CURRENT TYPE IS NOT IMPLEMENTED'
		STOP
	ENDIF ELSE BEGIN
		IF type EQ 'L' THEN dumarr = 'lonarr'
		IF type EQ 'D' THEN dumarr = 'dblarr'
		IF type EQ 'F' THEN dumarr = 'fltarr'
		IF type EQ 'L64' THEN dumarr = 'lon64arr'
	ENDELSE

	;;-----
	IF SIZE(input, /n_dimension) EQ 0L THEN input = [input]
	IF SIZE(input, /n_dimension) EQ 1L THEN n_dim = 1L
	IF SIZE(input, /n_dimension) GE 2L THEN n_dim = N_ELEMENTS(input(0,*))

	;;-----
	IF(n_dim EQ 1L) THEN n1 = N_ELEMENTS(input) 
	IF(n_dim GE 2L) THEN n1 = N_ELEMENTS(input(*,0)) 

	IF index LT 0L THEN BEGIN
		IF n_dim EQ 1L THEN BEGIN
			void	= EXECUTE('array = ' + STRTRIM(dumarr,2) + '(' + STRTRIM(unitsize,2) +  ')')
		ENDIF ELSE BEGIN
			void	= EXECUTE('array = ' + STRTRIM(dumarr,2) + '(' + $
				STRTRIM(unitsize,2) + ',' + STRTRIM(n_dim,2) + ')')
		ENDELSE
		nn = -1L
		n0	= 0L
		n1	= n1-1L
	ENDIF ELSE BEGIN
		nn = index
		n0	= nn
		n1	= nn + n1 -1L
	ENDELSE

	;;-----
	IF n_dim EQ 1L THEN n_size = N_ELEMENTS(array)
	IF n_dim GE 2L THEN n_size = N_ELEMENTS(array(*,0))

	;n0	= nn + 1L
	;n1	= nn + n1
	;n0	= nn
	;n1	= nn + n1 - 1L

	;;-----
	IF n1 GE n_size - 1L THEN BEGIN
		REPEAT BEGIN
			IF n_dim EQ 1L THEN void = EXECUTE('array = [array, ' + $
				STRTRIM(dumarr,2) + '(' + STRTRIM(unitsize,2)+ ')]')
	        	IF n_dim GE 2L THEN void = EXECUTE('array = [array, ' + $
				STRTRIM(dumarr,2) + '(' + STRTRIM(unitsize,2) + ',' + $
			        STRTRIM(n_dim,2) + ')]')
	
			IF n_dim EQ 1L THEN n_size = N_ELEMENTS(array)
			IF n_dim GE 2L THEN n_size = N_ELEMENTS(array(*,0))
		ENDREP UNTIL n1 LT n_size
	ENDIF

	IF n_dim EQ 1L THEN array(n0>0L:n1) = input
	IF n_dim GE 2L THEN array(n0>0L:n1,*) = input
	
	index = n1 + 1L
END


FUNCTION veluga::g_rotate, x, y, z, axis

	;;-----
	;; Rotate coordinate systems for axis to be z-axis
	;;
	;; If angular momentum is given, (x,y) is face-on and (x,z) or (y,z) are edge-on views
	;;-----

	rvec_z 	= axis / NORM(axis)
	const 	= SQRT( (1. - rvec_z(0)^2 - rvec_z(1)^2) / (1. - rvec_z(1)^2) )
	rvec_x 	= [const, 0., SQRT(1. - const^2)]
	rvec_y 	= CROSSP(rvec_z, rvec_x)

	xx0 	= x * rvec_x(0) + y * rvec_x(1) + z * rvec_x(2)
	yy0 	= x * rvec_y(0) + y * rvec_y(1) + z * rvec_y(2)
	zz0 	= x * rvec_z(0) + y * rvec_z(1) + z * rvec_z(2)
	;vx0 	= v * rvec_x(0) + v * rvec_x(1) + v * rvec_x(2)
	;vy0 	= v * rvec_y(0) + v * rvec_y(1) + v * rvec_y(2)
	;vz0 	= v * rvec_z(0) + v * rvec_z(1) + v * rvec_z(2)

	RETURN, {x:xx0, y:yy0, z:zz0}
END

FUNCTION veluga::g_boundind, xx=x, yy=y, zz=z, xr=xr, yr=yr, zr=zr
	tmp 	= 'ind = WHERE('

	IF KEYWORD_SET(xr) THEN BEGIN
		tmp += 'x GE xr(0) AND x LT xr(1) '
	ENDIF

	IF KEYWORD_SET(yr) THEN BEGIN
		tmp += 'AND y GE yr(0) AND y LT yr(1) '
	ENDIF

	IF KEYWORD_SET(zr) THEN BEGIN
		tmp += 'AND z GE zr(0) AND z LT zr(1) '
	ENDIF
	tmp     += ', nn)'
	void    = EXECUTE(tmp)

	RETURN, {ind:ind, n:nn}
END

;;----- Smoothing related
FUNCTION veluga::g_smooth_mafit, xx, yy2, nstep, dir, n_sigma
	yy 	= yy2
	nn 	= N_ELEMENTS(yy)

	CASE dir OF
		'F': BEGIN
			ind0 	= 0L
			ind1 	= nn-1L
			dn 		= 1L
			END
		'B': BEGIN
			ind0 	= nn-1L
			ind1 	= 0L
			dn 		= -1L
			END
		'N': BEGIN
			ind0 	= 0L
			ind1 	= nn-1L
			dn 	 	= 1L
			END
		ELSE: BEGIN
			self->errorout, 'wrong direction for mafit: check MA_direction'
			STOP
			END
	ENDCASE

    FOR i=ind0, ind1, dn DO BEGIN
    	i0 	= (i - LONG(nstep)/2) > 0L
		i1 	= (i + LONG(nstep)/2) < (nn-1L)

		IF i0 EQ 0L THEN i1 = (i + i - i0) < (nn-1L)
		IF i1 EQ nn-1L THEN i0 = (i - (i1-i)) > 0L

		IF dir EQ 'F' OR dir EQ 'B' THEN dummy 	= yy(i0:i1) ELSE dummy = yy2(i0:i1)
		


		IF ~KEYWORD_SET(sigma) THEN BEGIN
			yy(i) 	= MEAN(dummy)
		ENDIF ELSE BEGIN
			avg 	= MEAN(dummy)
			std 	= STDDEV(dummy)
			cut 	= WHERE( ABS(dummy - avg) LT std*n_sigma)
			yy(i) 	= MEAN( dummy(cut) )
		ENDELSE
	ENDFOR

	RETURN, yy
END

FUNCTION veluga::g_smooth, xx, yy, type=type, $
	MA_step=MA_step, MA_direction=MA_direction, MA_nsigma=MA_nsigma
	;;-----
	;; Line smoothing with different algorithm
	;;	xx, yy: [N] float/double
	;;
	;;  type: [1] string
	;;		'MA': moving average. MA_step & MA_direction should be argued
	;;
	;;
	;;	MA_step: [1] integer
	;;		filter width for computing average
	;;
	;;	MA_direction: [1] string
	;;		direction of filter
	;;		'F' or 'B' as forward or backward
	;;		'N' by no movements
	;;
	;;	MA_nsigma: [1] float/double
	;;		sigma-clipping for computing average
	;;		e.g., MA_nsigma=3 means that an average is computed with points within 3 sigma
	;;	
	;;
	;;-----

	IF ~KEYWORD_SET(type) THEN type = 'MA'
	type 	= STRUPCASE(type)

	CASE type OF
		'MA': BEGIN
			IF ~KEYWORD_SET(MA_step) THEN MA_step = 5L
			IF ~KEYWORD_SET(MA_direction) THEN MA_direction = 'F'
			IF ~KEYWORD_SET(MA_nsigma) THEN MA_nsigma = -1.d
			MA_direction 	= STRUPCASE(MA_direction)
			RETURN, self->g_smooth_mafit(xx, yy, MA_step, MA_direction, MA_nsigma)
			END
	ENDCASE

END
;;-----
;; SIMPLE GET FTNS
;;	-- RAMSES RELATED
;;-----
FUNCTION veluga::g_info, snap0

	;;-----
	;; Read info file
	;;	employed from rd_info.pro
	;;-----
	str 	= STRING(snap0, format='(I5.5)')
	settings 	= self.getheader()


  	my_narr		= LONARR(6)
	my_narr(0:5)  = 0L
	my_rarr       = DBLARR(11)
	my_rarr(0:10) = 0.
  	
	dataname      = string("",format='(a13)')

	file 	= settings.dir_raw + '/output_' + str + '/info_' + str + '.txt'
    OPENR,2,file
  
	value = 0L
	FOR il = 0,5 DO BEGIN
		READF, 2, dataname, value,format='(a13,I11)'
		my_narr(il) = value
	ENDFOR

	data="alors je sais vraiment pas qupi mettre pour que ca marche"

	READF, 2, data, format='(a50)'

	value =0d0
	FOR il = 0,10 DO BEGIN
		READF, 2, dataname, value, format='(a13,E23.15)'
		my_rarr(il) = value
	ENDFOR
	CLOSE, 2

	;; unit in cgs
	kpc     = 3.086d21;3.08568025e21
	twopi   = 6.2831853d0
	hplanck = 6.6262000d-27
	eV      = 1.6022000d-12
	kB      = 1.3806200d-16
	clight  = 2.9979250d+10
	Gyr     = 3.1536000d+16
	X       = 0.76
	Y       = 0.24 
	rhoc    = 1.8800000d-29
	mH      = 1.6600000d-24
	mu_mol  = 1.2195d0
	G       = 6.67259e-8
	m_sun   = 1.98892d33

	cgs 	= {kpc:kpc, hplanck:hplanck, eV:eV, kB:kB, clight:clight, Gyr:Gyr, mH:mH, G:G, m_sun:m_sun}

	scale_l    = my_rarr(8)
	scale_d    = my_rarr(9)
	scale_t    = my_rarr(10)

	;; scale_m convert mass in user units into g
	scale_m    = scale_d*(DOUBLE(scale_l))^3
	;; scale_v convert velocity in user units into cm/s
	scale_v    = scale_l / scale_t
	;; scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
	scale_T2   = mH/kB * scale_v^2.
	scale_nH   = X/mH * scale_d
	;; scale covert mettalicty into solar metallicity Z_sol = 0.02
	scale_Z    = 1./settings.sun_met 
	scale_flux = scale_v*scale_d*kpc*kpc*Gyr/m_sun
  
	info={cgs:cgs, boxtokpc:my_rarr(0)*scale_l/kpc,tGyr:my_rarr(1)*scale_t/Gyr,boxlen:my_rarr(0),levmin:my_narr(2),levmax:my_narr(3),unit_l:scale_l,unit_d:scale_d,unit_t:scale_t,unit_m:scale_m,unit_v:scale_v,unit_nH:scale_nH,unit_T2:scale_T2,unit_Z:scale_Z,kms:scale_v/1d5,unit_flux:scale_d*scale_v*(1e-9*Gyr)*(kpc)*(kpc)/m_sun,aexp:my_rarr(2), H0:my_rarr(3), omega_m:my_rarr(4), omega_l:my_rarr(5), omega_k:my_rarr(6), omega_b:my_rarr(7), ncpu:my_narr(0), ndim:my_narr(1)}
  
	;;-----
	;; Read Hilbert Indices
	;;-----
	ncpu	= my_narr(0)
	hindex	= DBLARR(ncpu,2)
	OPENR, 2, file
	str	= ' '
	FOR i=0L, 20L DO READF, 2, str
	FOR i=0L, ncpu-1L DO BEGIN
 	aa 	= DBLARR(3)
 	READF, 2, aa
 	hindex(i,0) = aa(1)
 	hindex(i,1) = aa(2)
	ENDFOR
	CLOSE, 2
	info	= CREATE_STRUCT(info,'hindex', hindex)

	;rd_info, info, file=settings.dir_raw + 'output_' + str + '/info_' + str + '.txt'
	RETURN, info
END

FUNCTION veluga::g_domain, snap0, xc2, yc2, zc2, rr2
	;;-----
	;; Return domain list
	;; 		positions are given in kpc unit
	;;-----

	settings	= self->getheader()
	info 	= self->g_info(snap0)

	n_gal 	= N_ELEMENTS(xc2)
	n_mpi 	= info.ncpu

	IF N_ELEMENTS(xc2) NE N_ELEMENTS(rr2) THEN rr2 = xc2*0.d + MAX(rr2)

	xc 	= DOUBLE(xc2)*info.cgs.kpc / info.unit_l
	yc 	= DOUBLE(yc2)*info.cgs.kpc / info.unit_l
	zc 	= DOUBLE(zc2)*info.cgs.kpc / info.unit_l
	rr 	= DOUBLE(rr2)*info.cgs.kpc / info.unit_l

	dom_list 	= LONARR(n_gal, n_mpi) - 1L

	ftr_name 	= settings.dir_lib + 'src/fortran/find_domain.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0) 	= n_gal
		larr(1) 	= n_mpi
		larr(2)		= 1L

		darr(0) 	= 1.d

	void 	= CALL_EXTERNAL(ftr_name, 'find_domain', $
		xc, yc, zc, rr, info.hindex, info.levmax, dom_list, larr, darr)

	void 	= WHERE(dom_list GE 0L, ncut)
	
	IF ncut GE 1L THEN BEGIN
		dom_all 	= LONARR(n_gal*n_mpi)-1L
		i0 	= 0L
		FOR i=0L, n_gal-1L DO BEGIN
			cut 		= WHERE(dom_list(i,*) GE 1L)
			cut2 		= (ARRAY_INDICES(dom_list(i,*), cut))(1,*)
			dom_list2	= self->g_unique(cut2) + 1L

			i1 	= i0 + N_ELEMENTS(dom_list2)-1L
			dom_all(i0:i1) 	= dom_list2
			i0 	= i1 + 1L
		ENDFOR
		dom_all 	= self->g_unique(dom_all)
		dom_all 	= dom_all(WHERE(dom_all GE 1L))
		RETURN, dom_all
	ENDIF ELSE BEGIN
		self->errorout, '		g_domain: NO DOMAIN LEFT'
		RETURN, 1L
	ENDELSE
END

FUNCTION veluga::g_part, snap0, xc2, yc2, zc2, rr2, dom_list=dom_list, simout=simout

	;;-----
	;; Read Particle within a sphere
	;; IF dom_list is set, read all ptcls in the argued domain list
	;;	/simout 	- output as the raw simulation unit
	;;-----

	settings	= self->getheader()
	num_thread = self.num_thread
	IF ~KEYWORD_SET(dom_list) THEN BEGIN
		dom_list 	= self->g_domain(snap0, xc2, yc2, zc2, rr2)
	ENDIF

	info 	= self->g_info(snap0)
	ncpu 	= N_ELEMENTS(dom_list)
	file 	= settings.dir_raw + 'output_' + STRING(snap0, format='(I5.5)') + '/part_' + STRING(snap0, format='(I5.5)') + '.out'
    TIC
	;;-----
	;; MEMORY ALLOCATION
	;;-----
	npart_tot 	= 0L
	part_ind 	= LONARR(ncpu)

	ftr_name 	= settings.dir_lib + 'src/fortran/jsrd_part_totnum.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0)	= ncpu
		larr(2) = num_thread
		larr(3)	= STRLEN(file)

	void 	= CALL_EXTERNAL(ftr_name, 'jsrd_part_totnum', $
		larr, darr, file, npart_tot, part_ind, dom_list)

    TOC, elapsed_time=elt1

    TIC
	xp 	= DBLARR(npart_tot,3)
	vp 	= DBLARR(npart_tot,3)
	mp 	= DBLARR(npart_tot)
	ap 	= DBLARR(npart_tot)
	zp 	= DBLARR(npart_tot)
	fam = LONARR(npart_tot)
	tag = LONARR(npart_tot)
	dl 	= DBLARR(npart_tot)
	id 	= LON64ARR(npart_tot)

	;;-----
	;; READ
	;;-----

	ftr_name 	= settings.dir_lib + 'src/fortran/jsrd_part.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0)	= ncpu
		larr(2)	= num_thread
		larr(3) = STRLEN(file)
		larr(4) = npart_tot

		;; skip domain reading
		larr(15) = 0L
		IF settings.skiprd_domain EQ 1L THEN larr(15) = 20L

		;; skip time reading
		larr(16) = 0L
		IF settings.skiprd_time EQ 1L THEN larr(16) = 20L

		;; skip metal reading
		larr(17) = 0L
		IF settings.skiprd_metal EQ 1L THEN larr(17) = 20L

		;; famtype ver
		larr(18) = 0L
		IF settings.famtype EQ 'new' THEN larr(18) = 20L

		;; id byte type
		larr(19) = 0L
		IF settings.idtype EQ 'long64' THEN larr(19) = 20L

	void 	= CALL_EXTERNAL(ftr_name, 'jsrd_part', $
		larr, darr, file, part_ind, xp, vp, mp, ap, zp, $
		fam, tag, dl, id, dom_list)

	part 	= self->allocate(npart_tot, type='part')
	;REPLICATE({xx:0.d, yy:0.d, zz:0.d, vx:0.d, vy:0.d, vz:0.d, mp:0.d, ap:0.d, zp:0.d, id:0L, family:0L, domain:0L}, npart_tot)

	IF ~KEYWORD_SET(simout) THEN BEGIN
		xp 	*= (info.unit_l/info.cgs.kpc)	;; [kpc]
		vp 	*= info.kms 					;; [kms]
		mp 	*= (info.unit_m / info.cgs.m_sun); [Msun]
	ENDIF

	part.xx 	= xp(*,0)
	part.yy 	= xp(*,1)
	part.zz 	= xp(*,2)

	part.vx 	= vp(*,0)
	part.vy 	= vp(*,1)
	part.vz 	= vp(*,2)

	part.mp 	= mp
	part.ap 	= ap
	part.zp 	= zp
	part.family = fam
	part.domain = dl
	part.id 	= id

	TOC, elapsed_time=elt2

    PRINT, elt1, ' - totnum ', elt2, ' - read all'
	RETURN, part
END

FUNCTION veluga::g_cell, snap0, xc2, yc2, zc2, rr2, dom_list=dom_list, simout=simout, timereport=timereport
	;;-----
	;; Read AMR cells within a sphere
	;;	snap0: [1] integer
	;;		snapshot number
	;; 
	;; 	xc2, yc2, zc2, rr2: [1] double
	;;		center & radius in kpc unit
	;;
	;;	dom_list: [N] integer
	;;		domain_list
	;;		If set, read all cells in the argued domain. If not, read all cells inside rr2 + dx
	;;
	;;	simout: boolean
	;;		If set, output as the raw code unit
	;;
	;;	timereport: boolean
	;;		If set, time report will be shown
	;;
	;;-----



	TIC
	settings	= self->getheader()
	num_thread = self.num_thread
	IF ~KEYWORD_SET(dom_list) THEN BEGIN
		dom_list 	= self->g_domain(snap0, xc2, yc2, zc2, rr2)
	ENDIF

	info 	= self->g_info(snap0)
	ncpu 	= N_ELEMENTS(dom_list)
	dir 	= settings.dir_raw + 'output_' + STRING(snap0, format='(I5.5)')

	xc	= xc2 / info.unit_l * info.cgs.kpc
	yc	= yc2 / info.unit_l * info.cgs.kpc
	zc	= zc2 / info.unit_l * info.cgs.kpc
	rr	= rr2 / info.unit_l * info.cgs.kpc

	IF ~KEYWORD_SET(xr) OR ~KEYWORD_SET(yr) OR ~KEYWORD_SET(zr) THEN BEGIN
		xr	= [-1d,1d] * rr + xc
		yr	= [-1d,1d] * rr + yc
		zr	= [-1d,1d] * rr + zc
	ENDIF

	IF ~KEYWORD_SET(range_refine) THEN range_refine = [-1.d, 2.d]
	IF ~KEYWORD_SET(range_xx) THEN range_xx = [-2.d, 2.d]
	IF ~KEYWORD_SET(range_yy) THEN range_yy = [-2.d, 2.d]
	IF ~KEYWORD_SET(range_zz) THEN range_zz = [-2.d, 2.d]

	file_a	= dir + '/amr_' + STRING(snap0,'(I5.5)') + '.out'
	file_h	= dir + '/hydro_' + STRING(snap0,'(I5.5)') + '.out'
	file_i	= dir + '/info_' + STRING(snap0,'(I5.5)') + '.txt'

	TOC, elapsed_time=time_header
	TIC
	;;-----
	;; MEMORY ALLOCATE
	;;-----
	ftr_name 	= settings.dir_lib + 'src/fortran/jsamr2cell_totnum.so'
	larr = LONARR(20) & darr = DBLARR(20)
		larr(0) = ncpu
		larr(2) = 1L;settings.num_thread
		larr(3) = STRLEN(file_a)
		larr(4) = STRLEN(file_h)
		larr(5) = LONG(info.ncpu)
		larr(6) = LONG(info.ndim)
		larr(7) = LONG(info.levmin)
		larr(8) = LONG(info.levmax)

		mg_ind	= LONARR(info.ncpu)
		ntot	= 0L
		nvarh	= 0L

		
		void	= CALL_EXTERNAL(ftr_name, 'jsamr2cell_totnum', $
			larr, darr, file_a, file_h, ntot, nvarh, mg_ind, dom_list)

	IF ntot EQ 0L THEN BEGIN
		self->errorout,'No leaf cells in this domain. Dummy array returned'
		cell 	= self->allocate(1L, type='cell')
		cell.xx 	= -1.d
		RETURN, cell
	ENDIF

	mesh_xg	= DBLARR(ntot,info.ndim)
	mesh_vx	= DBLARR(ntot,info.ndim)
	mesh_dx	= DBLARR(ntot)
	mesh_hd	= DBLARR(ntot,nvarh)
	mesh_lv	= LONARR(ntot)*0L - 10L
	mesh_mp	= DBLARR(ntot)

	nx 	= larr(10)
	ny 	= larr(11)
	nz 	= larr(12)
	nboundary 	= larr(13)

	levelind= LONARR(info.levmax,3) ;; for memory efficiency


	;;-----
	;; READ CELL
	;;-----
	ftr_name	= settings.dir_lib + 'src/fortran/jsamr2cell.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0)	= ncpu;icpu
		larr(2)	= self.num_thread
		larr(3)	= STRLEN(file_a)
		larr(4)	= STRLEN(file_h)
		larr(5)	= STRLEN(file_i)
		larr(6) = LONG(info.ncpu)
		larr(7) = LONG(info.ndim)
		larr(8) = LONG(info.levmin)
		larr(9) = LONG(info.levmax)
		larr(10)= ntot
		larr(11)= nvarh
		larr(12)= nx
		larr(13)= ny
		larr(14)= nz
		larr(15)= nboundary

		IF ~KEYWORD_SET(simout) THEN BEGIN
			tokpc 	= (info.unit_l/info.cgs.kpc)	;; [kpc]
			tokms	= info.kms 						;; [km/s]
			tomsun 	= info.unit_d * info.unit_l^3.d / info.cgs.m_sun 	;; [Msun]
			tocc 	= info.unit_nH 					;; [/cc]
			toKmu	= info.unit_T2 					;; [K/mu]
		ENDIF ELSE BEGIN
			tokpc 	= 1.d
			tokms 	= 1.d
			tomsun 	= 1.d
			tocc 	= 1.d
			toKmu 	= 1.d
		ENDELSE

		darr(0)	= tokpc
		darr(1)	= tokms
		darr(2)	= tomsun
		darr(3)	= tocc
		darr(4)	= toKmu

		void	= CALL_EXTERNAL(ftr_name, 'jsamr2cell', $
			larr, darr, file_a, file_h, file_i, $
			mg_ind, mesh_xg, mesh_dx, mesh_hd, mesh_lv, mesh_mp, dom_list, levelind)


	TOC, elapsed_time=time_read

	TIC
	;;-----
	;; POST PROCESSING
	;;-----
	;cut	= WHERE(mesh_lv GE 0L,ncell)
	;IF ncell NE N_ELEMENTS(mesh_lv) THEN BEGIN
	;	self->errorout, '!?'
	;	STOP
	;ENDIF
	cell 	= self->allocate(ntot, type='cell')

	TOC, elapsed_time=time_allocation
	TIC
	;lind 	= REPLICATE({i0:0L, i1:0L}, info.levmax+1L)
	;lind(1:*).i0 	= levelind(*,0)
	;lind(1:*).i1 	= levelind(*,1)

	cell.xx 	= mesh_xg(*,0)
	cell.yy 	= mesh_xg(*,1)
	cell.zz 	= mesh_xg(*,2)

	cell.vx 	= mesh_hd(*,1)
	cell.vy 	= mesh_hd(*,2)
	cell.vz 	= mesh_hd(*,3)

	cell.level 	= mesh_lv
	cell.dx 	= mesh_dx

	cell.den 	= mesh_hd(*,0)
	cell.temp 	= mesh_hd(*,4)
	cell.zp 	= mesh_hd(*,5)
	cell.mp 	= mesh_mp

	cell.UE 	= mesh_hd(*,4)/(5.d/3.-1.d) * info.unit_T2 / (1.66d-24) * 1.38049d-23 * 1e-3 / toKmu ;; [km/s]^2
	cell.p_thermal = mesh_hd(*,4)*mesh_hd(*,0) * info.unit_m / info.unit_l / info.unit_t^2 / 1.3806200d-16

	cell.levelind 	= [levelind(0,*), levelind]

	IF N_ELEMENTS(settings.hydro_variables) GT 7L THEN BEGIN
		FOR i2=6L, N_ELEMENTS(settings.hydro_variables)-1L DO BEGIN
			IF STRPOS(settings.hydro_variables(i2),'skip') GE 0L THEN CONTINUE

			str 	= 'cell.' + STRTRIM(settings.hydro_variables(i2),2) + $
				' = mesh_hd(*,' + STRING(i2) + ')'
			void 	= EXECUTE(str)
			;; [mass frac]
		ENDFOR
	ENDIF
	
	TOC, elapsed_time=time_post

;	mesh_xg 	= mesh_xg(cut,*)
;	mesh_hd 	= mesh_hd(cut,*)
;	mesh_lv 	= mesh_lv(cut)
;	mesh_dx 	= mesh_dx(cut)
;
;	cell.xx 	= mesh_xg(*,0)
;	cell.yy 	= mesh_xg(*,1)
;	cell.zz 	= mesh_xg(*,2)
;
;	cell.vx 	= mesh_hd(*,1)
;	cell.vy 	= mesh_hd(*,2)
;	cell.vz 	= mesh_hd(*,3)
;
;	cell.level 	= mesh_lv
;	cell.dx 	= mesh_dx
;
;	cell.den 	= mesh_hd(*,0)
;	cell.temp 	= mesh_hd(*,4)
;	cell.zp 	= mesh_hd(*,5)
;	;cell.levelind 	= PTR_NEW(levelind)
;
;	IF N_ELEMENTS(settings.hydro_variables) GT 7L THEN BEGIN
;		FOR i=6L, N_ELEMENTS(settings.hydro_variables)-1L DO BEGIN
;			IF STRPOS(settings.hydro_variables(i),'skip') GE 0L THEN CONTINUE
;
;			str 	= 'cell.' + STRTRIM(settings.hydro_variables(i),2) + ' = mesh_hd(*,' + STRING(i) + ')'
;			void 	= EXECUTE(str)
;			;; [mass frac]
;		ENDFOR
;	ENDIF
;
;	cell.UE 	= mesh_hd(*,4)/mesh_hd(*,0)/(5.d/3.-1.d) * info.unit_T2 / (1.66d-24) * 1.38049d-23 * 1e-3 ;; [km/s]^2
;	cell.p_thermal = mesh_hd(*,4) * info.unit_m / info.unit_l / info.unit_t^2 / 1.3806200d-16
	IF KEYWORD_SET(timereport) THEN BEGIN
		PRINT, '%%----- TIME REPORT FOR G_CELL'
		PRINT, '	READ HEADER : ', time_header, '[sec]'
		PRINT, '	READ CELL 	: ', time_read, '[sec]'
		PRINT, '	ALLOCATION 	: ', time_allocation, '[sec]'
		PRINT, '	COPY ARRAY 	: ', time_post, '[sec]'
	ENDIF
	RETURN, cell

END


FUNCTION veluga::g_cfrac, snap0, xc, yc, zc, aperture
	;;-----
	;; Compute contamination fractions
	;;		positions are given in kpc unit
	;;		aperture in kpc unit
	;;			shuld be in a [N_gal, N_aper] format
	;;-----

	settings	= self->getheader()
	num_thread 	= self.num_thread
	
	n_gal 	= N_ELEMENTS(xc)
	n_aper	= N_ELEMENTS(aperture(0,*))

	;;------
	;; Get Domain
	;;------
	rr 	= xc
	FOR i=0L, n_gal-1L DO rr(i) = MAX(aperture(i,*))

	dom_all 	= self->g_domain(snap0, xc, yc, zc, rr)


	;;-----
	;; Read all ptcls
	;;-----
	part 	= self->g_part(snap0, 0.d, 0.d, 0.d, 0.d, dom_list=dom_all)

	dm_ind 	= WHERE(part.family EQ 1L, nn_dm)
	part 	= self->g_extract(part, dm_ind)
	xp 	= DBLARR(nn_dm,3)
	xp(*,0)	= part.xx
	xp(*,1)	= part.yy
	xp(*,2)	= part.zz
	mp 		= part.mp

	;;-----
	;; CFrac computation
	;;-----
	info 	= self->g_info(snap0)

	conf_n 	= DBLARR(n_gal, n_aper)
	conf_m 	= DBLARR(n_gal, n_aper)

	xc0 	= xc; * 3.086d21/info.unit_l
	yc0 	= yc; * 3.086d21/info.unit_l
	zc0 	= zc; * 3.086d21/info.unit_l
	aperture= aperture; * 3.086d21/info.unit_l
	;rr0 	= rr * 3.086d21/info.unit_l

	dmp_mass 	= 1./(settings.neff*1.d)^3 * (info.omega_M - info.omega_b) / info.omega_m
	dmp_mass 	*= (info.unit_m / info.cgs.m_sun)

	ftr_name 	= settings.dir_lib + 'src/fortran/get_contam.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0) = n_gal
		larr(1) = nn_dm
		larr(2) = n_aper
		larr(3) = num_thread
		larr(4) = 1L
		larr(5)	= 64L
		larr(11)= STRLEN(settings.dir_raw)


		darr(0)	= dmp_mass

	void 	= CALL_EXTERNAL(ftr_name, 'get_contam', $
		larr, darr, settings.dir_raw, xc0, yc0, zc0, aperture, $
		xp, mp, conf_n, conf_m)

	RETURN, {N:conf_n, M:conf_m}
END

FUNCTION veluga::g_gyr, snap0, tconf
	;;-----
	;; Compute age of ptcls in [Gyr, redsh, scale factor] from conformal unit at a given snapshot
	;;-----
	settings= self->getheader()
	info 	= self->g_info(snap0)
	oM	= info.omega_m
	oL	= info.omega_l
	H0	= info.h0

	;;-----
	;; Run or Load - conformal time table
	;;-----

	conf_tbl 	= self->t_conformal_load(oM, oL)
	lbt_tbl 	= self->t_lbt_load(oM, oL, H0)

	;;-----
	;; INTERPOLATION
	;;-----
	n_part	= N_ELEMENTS(tconf)
	;t_res 	= DBLARR(n_part, 2) - 1.0d8 ;; [SFactor, GYR]
	v1 	= DBLARR(n_part) - 1.0d8 
	v2 	= DBLARR(n_part) - 1.0d8

	ftr_name	= settings.dir_lib + 'src/fortran/prop_time.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0)	= n_part
		larr(1)	= self.num_thread

		darr(0)	= 1./info.aexp - 1.0d

	void	= CALL_EXTERNAL(ftr_name, 'prop_time', $
		tconf, v1, v2, conf_tbl.conft, conf_tbl.sfact, lbt_tbl.redsh, lbt_tbl.gyr, $
		larr, darr)

	age 	= REPLICATE({gyr:0.d, redsh:0.d, sfact:0.d}, n_part)
	age.sfact 	= v1
	age.redsh 	= 1./v1 - 1.d
	age.gyr 	= v2

	RETURN, age
END

FUNCTION veluga::g_ztoage, z, info=info

	IF ~KEYWORD_SET(info) THEN info = self->g_info(1L)
	lbt_tbl 	= self->t_lbt_load(info.omega_m, info.omega_l, info.h0)


	RETURN, lbt_tbl.gyr(-1) - INTERPOL(lbt_tbl.gyr, lbt_tbl.redsh, z)
END
FUNCTION veluga::g_sfr, xx, yy, zz, age, mass, xc, yc, zc, aperture=aperture, timewindow=timewindow
	;;-----
	;; Compute SFR with given ptcls
	;;	- xx, yy, and zz
	;;		positions of ptcls in kpc unit (physical)
	;;	- age, mass
	;;		age of ptcls in Gyr
	;;		mass of ptcls in Msun
	;;	- xc, yc, and zc
	;;		positions of the center
	;;	- aperture
	;;		aperture size in kpc unit
	;;		if set to be negative, use all ptcls
	;;	- timewindow
	;;		timewindow for SFR mesurement in Gyr unit
	;;-----


	;; Default
	IF ~KEYWORD_SET(aperture) THEN apreture = 10.d
	IF ~KEYWORD_SET(timewindow) THEN timewindow = 0.1d
	settings= self->getheader()

	;; Allocate
	n_part 	= N_ELEMENTS(xx)

	ftr_name	= settings.dir_lib + 'src/fortran/prop_sfr.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0)	= 1L
		larr(1)	= n_part
		larr(2)	= self.num_thread

		darr(0)	= aperture
		darr(10)= timewindow

		void	= CALL_EXTERNAL(ftr_name, 'prop_sfr', $
			larr, darr, $
			DOUBLE(xx), DOUBLE(yy), DOUBLE(zz), DOUBLE(age), DOUBLE(mass), $
			DOUBLE(xc), DOUBLE(yc), DOUBLE(zc))

	sfr = darr(19)
	RETURN, sfr
END

FUNCTION veluga::g_luminosity, mp, ap, zp, band

	;;-----
	;; Get Luminosity of given ptcls
	;;
	;;	mp in Msun
	;; 	ap in Gyr
	;;  zp in Zsun
	;;
	;;	result, Luminosity of ptcls in Lsun
	;;-----
	tbl_sdss 	= self->t_miles_sdss_load()
	tbl_galex 	= self->t_miles_galex_load()

	flux 	= REPLICATE({u:0.d, g:0.d, r:0.d, i:0.d, z:0.d, nuv:0.d}, N_ELEMENTS(mp))

	FOR i=0L, N_ELEMENTS(band)-1L DO BEGIN
		;; Set Table
		CASE band(i) OF
			'u'		: BEGIN
				ref_ml = tbl_sdss.u & ref_met = tbl_sdss.metal & ref_age = tbl_sdss.age & mag0 = 6.55 & END
			'g'		: BEGIN
				ref_ml = tbl_sdss.g & ref_met = tbl_sdss.metal & ref_age = tbl_sdss.age & mag0 = 5.12 & END
			'r'		: BEGIN
				ref_ml = tbl_sdss.r & ref_met = tbl_sdss.metal & ref_age = tbl_sdss.age & mag0 = 4.68 & END
			'i'		: BEGIN
				ref_ml = tbl_sdss.i & ref_met = tbl_sdss.metal & ref_age = tbl_sdss.age & mag0 = 4.57 & END
			'z'		: BEGIN
				ref_ml = tbl_sdss.z & ref_met = tbl_sdss.metal & ref_age = tbl_sdss.age & mag0 = 4.54 & END
			'nuv'	: BEGIN
				ref_ml = tbl_galex.nuv & ref_met = tbl_galex.metal & ref_age = tbl_galex.age & mag0 = 10.18 & END
			'NUV'	: BEGIN
				ref_ml = tbl_galex.nuv & ref_met = tbl_galex.metal & ref_age = tbl_galex.age & mag0 = 10.18 & END
			ELSE: STOP
		ENDCASE

		;; Extrapolation & Rearrange
		age_arr = ap & met_arr = zp

		age_arr = age_arr > MIN(ref_age) & age_arr = age_arr < MAX(ref_age)
		met_arr = met_arr > MIN(met_arr) & met_arr = met_arr < MAX(met_arr)

			
		;; Interpolation
		age_ind = LONARR(N_ELEMENTS(age_arr))-1L
		met_ind = LONARR(N_ELEMENTS(age_arr))-1L

		FOR j=1L, N_ELEMENTS(ref_age)-1L DO BEGIN
			cut 	= WHERE(age_ind EQ -1L AND age_arr - ref_age(j) LT 0., ncut)
			IF ncut GE 1L THEN age_ind(cut)	= j-1
		ENDFOR
		cut 	= WHERE(age_arr GT MAX(ref_age), ncut)
		IF ncut GE 1L THEN age_ind(cut) = N_ELEMENTS(ref_age)-2L

		FOR j=1L, N_ELEMENTS(ref_met)-1L DO BEGIN
			cut 	= WHERE(met_ind EQ -1L AND met_arr - ref_met(j) LT 0., ncut)
			IF ncut GE 1L THEN met_ind(cut)	= j-1
		ENDFOR
		cut 	= WHERE(met_arr GT MAX(ref_met), ncut)
		IF ncut GE 1L THEN met_ind(cut) = N_ELEMENTS(ref_met)-2L

		z00 	= ref_ml(age_ind, met_ind)
		z01 	= ref_ml(age_ind, met_ind+1L)
		z10 	= ref_ml(age_ind+1L, met_ind)
		z11 	= ref_ml(age_ind+1L, met_ind+1L)

		zz0 	= (z10-z00) / (ref_age(age_ind+1L) - ref_age(age_ind)) * (age_arr - ref_age(age_ind)) + z00
		zz1 	= (z11-z01) / (ref_age(age_ind+1L) - ref_age(age_ind)) * (age_arr - ref_age(age_ind)) + z01

		ml 		= (zz1-zz0) / (ref_met(met_ind+1L) - ref_met(met_ind)) * (met_arr - ref_met(met_ind)) + zz0

		lu 	= 1.d/ml * mp ;; in Lsun



		CASE band(i) OF
			'u'		: flux.u 	= lu
			'g'		: flux.g 	= lu
			'r'		: flux.r 	= lu
			'i'		: flux.i 	= lu
			'z'		: flux.z 	= lu
			'nuv'	: flux.nuv 	= lu
			'NUV'	: flux.nuv 	= lu
		ENDCASE
	ENDFOR

	RETURN, flux
END

FUNCTION veluga::g_mag, lum, band

	;;-----
	;; based on the MILES catalog with assuming that L/Lsun given in MILES uses Lsun in each band not a bolometric value
	;;-----

	FOR i=0L, N_ELEMENTS(band)-1L DO BEGIN
		CASE band(i) OF
			'u'		: magsun = 6.55
			'g'		: magsun = 5.12
			'r'		: magsun = 4.68
			'i'		: magsun = 4.57
			'z'		: magsun = 4.54
			'nuv'	: magsun = 10.18
			'NUV'	: magsun = 10.18
		ENDCASE

		mag 	= magsun - 2.5 * ALOG10(TOTAL(lum))
	ENDFOR

	RETURN, mag
END

FUNCTION veluga::g_sbf, lum, size, band

	;;-----
	;; L_sun in each band = L0 * 10^(-M / 2.5) where L0 is 3.0128e28
	;;
	;;	lum in Lsun
	;;	size in kpc
	;;-----

	FOR i=0L, N_ELEMENTS(band)-1L DO BEGIN
		CASE band(i) OF
			'u'		: magsun = 6.55
			'g'		: magsun = 5.12
			'r'		: magsun = 4.68
			'i'		: magsun = 4.57
			'z'		: magsun = 4.54
			'nuv'	: magsun = 10.18
			'NUV'	: magsun = 10.18
		ENDCASE

		sbf_in_physical 	= TOTAL(lum) / (size*1d3)^2 	;; [Lsun / pc^2]
		sbf 	= magsun + 21.572 - 2.5d * ALOG10(sbf_in_physical)
	ENDFOR

	RETURN, sbf
END

FUNCTION veluga::g_potential, xx, yy, zz, mm, $
	p_type=p_type, e_type=e_type, bsize=bsize

	;;-----
	;; Compute potential using all mass components
	;;
	;;		PE in [(km/s)^2]
	;;-----

	settings 	= self->getheader()
	IF ~KEYWORD_SET(p_type) THEN p_type = 'mesh'
	IF ~KEYWORD_SET(e_type) THEN e_type = 'pole'
	IF ~KEYWORD_SET(bsize) THEN bsize = 128L

	;;-----
	;; READ PARTICLES
	;;-----
	
	;d3d_part 	= self->g_d3d(part.xx, part.yy, part.zz, [x0, y0, z0])
	;cut_part 	= WHERE( (part.family EQ 1L OR part.family EQ 2L) AND d3d_part LT r0, npart)

	;d3d_cell 	= self->g_d3d(cell.xx, cell.yy, cell.zz, [x0, y0, z0])
	;cut_cell 	= WHERE(d3d_cell LT r0, ncell)

	;;-----
	;; ALLOCATE
	;;-----
	npart 	= N_ELEMENTS(xx)
	pos 	= DBLARR(npart,3)
	mass 	= DOUBLE(mm)

	pos(*,0) 	= xx
	pos(*,1) 	= yy
	pos(*,2) 	= zz

	pot 	= DBLARR(npart)
	force 	= DBLARR(npart)
	IF N_ELEMENTS(mass) LE bsize THEN bsize = 4L

	;;-----
	;; Compute potential
	;;-----
	Gconst 		= 6.67408d-11 		;; m^3 kg^-1 s^-2
	mtokpc 		= (1./ 3.086d19)
	kgtoMsun 	= (1./1.98892d30)
	Gconst 		*= (mtoKpc / kgtoMsun * 1d-6)	;; (km/s)^2 Kpc/Msun


	ftr_name 	= settings.dir_lib + 'src/fortran/js_getpt_ft.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0)	= N_ELEMENTS(mass)
		larr(1) = 3L	;; dimension
		larr(2)	= self.num_thread
		IF p_type EQ 'mesh' THEN larr(3) = 0L ELSE IF p_type EQ 'pm' THEN larr(3) = 1L
		IF e_type EQ 'pole' THEN larr(4) = 0L ELSE IF e_type EQ 'part' THEN larr(4) = 1L

		larr(10)	= 0L 	;; Tree dimension spliting type
		larr(11)	= 0L 	;; Tree value spliting type
		larr(12) 	= bsize

		darr(0)	= Gconst

	void 	= CALL_EXTERNAL(ftr_name, 'js_getpt_ft', $
		larr, darr, pos, mass, pot, force)

	;part.PE 	= pot(0L:npart-1L) 		;; [km/s]^2
	;cell.PE 	= pot(npart:npart+ncell-1L)

	;part(cut_part).KE 	= self->g_d3d( part(cut_part).vx, part(cut_part).vy, part(cut_part).vz, []
	RETURN, {PE:pot, force:force}
END

FUNCTION veluga::g_extract, array, ind
	;;-----
	;; Reshape a particle / cell array with the argued ind
	;;
	;; For a cell array, levind is changed with the new shape
	;;-----

	tag 	= TAG_NAMES(array)

	n_old 	= array.N
	n_new 	= N_ELEMENTS(ind)

	IF array.tag EQ 'cell' THEN BEGIN
		levind 	= array.levelind * 0L

		cell_lev 	= array.level(ind)
		info 	= self->g_info(1L)

		ind0	= 0L
		FOR lev=0L, info.levmax DO BEGIN
			cut 	= WHERE(cell_lev EQ lev, ncut)


			IF ncut GE 1L THEN BEGIN
				ind1 	= ind0 + ncut - 1L
				levind(lev,0) 	= ind0
				levind(lev,1)	= ind1
				levind(lev,2) 	= ncut

				ind0 	= ind1 + 1L
			ENDIF
		ENDFOR
	ENDIF

	strdum 	= 'array2 = {'
	FOR i=0L, N_ELEMENTS(tag)-1L DO BEGIN
		CASE tag(i) OF
			'N': strdum = strdum + 'N:' + STRTRIM(n_new,2) + ' '
			'LEVELIND': strdum = strdum + 'levelind:levind '
			'TAG': strdum = strdum + 'tag:array.tag '
			ELSE: strdum = strdum + tag(i) + ':array.' + tag(i) + '(ind) '
			
		ENDCASE

		IF i NE N_ELEMENTS(tag)-1L THEN strdum = strdum + ', '
	ENDFOR
	strdum 	= strdum + '}'

	void 	= EXECUTE(strdum)

	RETURN, array2

END

FUNCTION veluga::g_celltype, n_snap, cell, xc, yc, zc, rc, vxc, vyc, vzc, dom_list=dom_list, n_shell=n_shell, bsize=bsize

	;;-----
	;; Get Cell type (ISM, CGM, IGM) around the given center
	;; 		Cell is classified based on Rhee+24
	;; 		PE & KE are measured and stored in the cell array
	;;		Particles (< rc) are used to compute potential
	;;		Cells outside the paerture are classified as IGM (and are not used for potential calculation)
	;; 
	;;
	;; Result: [N] integer
	;;		N is the number of cell
	;;		1  ISM
	;;		0  CGM
	;; 		-1 IGM (or surrounding cell)
	;;	
	;;	n_snap: [1] integer
	;;		snapshot number
	;;
	;;	cell: [] cell array
	;;
	;; 	xc, yc, zc, rr: [1] double
	;;		center & radius in kpc unit
	;;
	;;	vxc, vyc, vzc: [1] double
	;;		velocity of the center to get the relative velocities
	;;
	;;	dom_list: [N] integer
	;;		domain_list
	;;		If set, read all cells in the argued domain. If not, read all cells inside rr2 + dx
	;;
	;;	n_shell: [1] integer
	;;		# of radial bins when computing metallicity radial distribution of ISM to distinguish IGM / CGM
	;;
	;;	bsize: [1] integer
	;;		Leaf node size of KDTree when computing potential
	;;-----

	IF ~KEYWORD_SET(dom_list) THEN BEGIN
		dom_list 	= self->g_domain(n_snap, xc, yc, zc, rc*2.d)
	ENDIF

	IF ~KEYWORD_SET(n_shell) THEN n_shell = 100L
	IF ~KEYWORD_SET(bsize) THEN bsize = 1024L

	c_d3d 	= self->g_d3d(cell.xx, cell.yy, cell.zz, [xc, yc, zc])
	c_outside 	= WHERE(c_d3d GT rc, nc_outside)

	;;----- Read Part	
	part 	= self->g_part(n_snap, 0.d, 0.d, 0.d, 0.d, dom_list=dom_list)

	;;----- Get Part within the aperture
	p_d3d 	= self->g_d3d(part.xx, part.yy, part.zz, [xc, yc, zc])
	p_ind	= WHERE( (part.family EQ 1L OR part.family EQ 2L) AND p_d3d LT rc , np)
	IF np EQ 0L THEN STOP
	part 	= self->g_extract(part, p_ind)


	
	;;-----
	;; Get Potential
	;;-----
	nc 	= cell.n
	nn 	= np + nc
	dumx 	= DBLARR(nn)
	dumy 	= DBLARR(nn)
	dumz 	= DBLARR(nn)
	dumm 	= DBLARR(nn)

	dumx(0L:nc-1L)	 	= cell.xx
	dumx(nc:np+nc-1L)	= part.xx

	dumy(0L:nc-1L)	 	= cell.yy
	dumy(nc:np+nc-1L)	= part.yy

	dumz(0L:nc-1L)	 	= cell.zz
	dumz(nc:np+nc-1L)	= part.zz

	dumm(0L:nc-1L)	 	= cell.mp
	IF nc_outside GE 1L THEN dumm(c_outside) = 0.d

	dumm(nc:np+nc-1L)	= part.mp

	pot 	= self->g_potential(dumx, dumy, dumz, dumm, bsize=bsize)
	
	cell.PE 	= pot.PE(0L:nc-1L)
	;IF nc_outside GE 1L THEN cell.PE(c_outside) = 0.d

	cell.KE 	= 0.5d * (self->g_d3d(cell.vx, cell.vy, cell.vz, [vxc, vyc, vzc]))^2

	Etot 	= cell.PE + cell.KE + cell.UE
	
	;;-----
	;; Cell type with metallicity condition
	;;	1 : ISM
	;;	0 : CGM
	;;	-1: IGM
	;;-----
	cell_type 	= LONARR(N_ELEMENTS(Etot)) - 1L

	Ecut 	= 0.d 		;; Energy cut for boundness
	d_shell	= rc/n_shell
	minZval	= 0.1*0.02 	;; Lower bound Metallicity for CGM (0.1 Zsun)
	tempcut	= 1e7

	;; ISM by Bound cells
	cut 	= WHERE(Etot LT Ecut, ncut)
	IF ncut GE 1L THEN cell_type(cut) = 1L

	;; Compute radial metallicity distribution of ISM
	ism_met	= DBLARR(n_shell,2)
		ism_met(*,0)	= 1e8
		ism_met(*,1)	= 0.


	FOR i=0L, n_shell-1L DO BEGIN
		r0 	= d_shell * i
		r1 	= d_shell * (i+1.d)
		
		cut = WHERE(c_d3d GE r0 and c_d3d LT r1 AND cell_type EQ 1L, ncut)
		IF ncut EQ 0L THEN CONTINUE

		ism_met(i,0)	= self->g_wmean(cell.zp(cut), cell.mp(cut))
		ism_met(i,1) 	= self->g_wstddev(cell.zp(cut), cell.mp(cut))
	ENDFOR

	;; Extrapolate ism_met beyond the ISM boundary
	;;	, to keep high metallicity
	FOR i=1L, n_shell-1L DO BEGIN
		IF ism_met(i,0) GE 1e7 THEN BEGIN
			ism_met(i,*)	= ism_met(i-1,*)
		ENDIF
	ENDFOR

	;; CGM by 1) Positive E & 2) Z > MAX(Z_ism - dz_ism, minZval) 3) Outflowing or T > 1e7K
	vdot 	= (cell.xx - xc) * (cell.vx - vxc) + (cell.yy - yc) * (cell.vy - vyc) + (cell.zz - zc) * (cell.vz - vzc)
	vdot0	= MAX(ABS(vdot)) * (-2.d)

	FOR i=0L, n_shell-1L DO BEGIN
		r0 	= d_shell*i
		r1 	= d_shell*(i+1.d)

		cut = WHERE(c_d3d GE r0 AND c_d3d LT r1 AND cell_type NE 1L, ncut) ;; here 1) is already satisfied
		IF ncut EQ 0L THEN CONTINUE


		met_avg	= self->g_wmean(cell.zp(cut), cell.mp(cut))
		met_std = self->g_wstddev(cell.zp(cut), cell.mp(cut))

		lowZval 	= MAX([minZval, ism_met(i,0) - ism_met(i,1)])


		cut2 	= WHERE($
			(c_d3d GE r0 AND c_d3d LT r1) AND $
			cell_type NE 1L AND $ 		;; Condition 1)
			cell.zp GT lowZval AND $	;; Condition 2)
			(vdot GT 0. OR cell.temp GT tempcut) $ ;; Condition 3)
			, nc2)

		IF nc2 GE 1L THEN cell_type(cut2) = 0L
	ENDFOR

	;; IGM for left cells
	cell_type(c_outside) 	= -1L
	
	RETURN, cell_type
END

FUNCTION veluga::g_cellphase, n_snap, cell

	;;-----
	;; Get Cell phase based on cold (warm) vs. hot based on the criterion by Torrey+12
	;;
	;; Result: [N] integer
	;;		N is the number of cell
	;;		2  SF (nH > 10 cc)
	;;		1  cold
	;;		-1 hot
	;;	
	;;	n_snap: [1] integer
	;;		snapshot number
	;;
	;;	cell: [] cell array
	;;
	;;		density and temperature should be given in K and cc unit, respectively, which is the default output of g_cell
	;;-----

	info 	= self->g_info(n_snap)

	den2    = ALOG10(cell.den) + ALOG10(1.6600000d-24)   ;; g/cc
        den2    -= ALOG10(info.cgs.m_sun)   ;; Msun/cc
        den2    += 3.*ALOG10(info.cgs.kpc)  ;; Msun/Kpc^3
        den2    -= 2.*ALOG10(info.H0/100.)     ;; Msun h^2 / Kpc^3
        den2    -= ALOG10(1e10)         ;; 1e10Msun h^2 / Kpc^3
        den2    = 10.d^den2

    phase	= LONARR(cell.n) - 1L
    cut_cold= WHERE(ALOG10(cell.temp) LT 6. + 0.25 * ALOG10(den2), nc)
    IF nc GE 1L THEN phase(cut_cold) = 1L

    cut_sf 	= WHERE(cell.den GT 10., nsf)
    IF nsf GE 1L THEN phase(cut_sf) = 2L

    RETURN, phase
END

;;-----
;; DRAWING ROUTINES
;;-----
FUNCTION veluga::d_minmax, map2, min2, max2, stype=stype, loga=loga
	map 	= map2

	IF ~KEYWORD_SET(stype) THEN stype = 'log'
	IF ~KEYWORD_SET(loga) THEN loga = 1000.d

	CASE stype OF
		'log' : BEGIN
			map 	= ALOG(loga*map + 1.d) / ALOG(loga)
			min 	= ALOG(loga*min2 + 1.d) / ALOG(loga)
			max 	= ALOG(loga*max2 + 1.d) / ALOG(loga)
			END
		'log2' : BEGIN
			map 	= ALOG10(map + 1.d)
			min 	= ALOG10(min2 + 1.d)
			max 	= ALOG10(max2 + 1.d)
			END
		'lin' : BEGIN
			min = min2 & max = max2

			END
	ENDCASE

	cut 	= WHERE(map LT min, ncut)
	IF ncut GE 1L THEN map(cut) = 0.

	cut 	= WHERE(map GT 0.)
	map(cut)-= min

	cut 	= WHERE(map GT max-min, ncut)
	IF ncut GE 1L THEN map(cut) = max-min

	map 	= map / (max - min)

	RETURN, map

END


FUNCTION veluga::d_2dmap, xx, yy, zz=zz, xr=xr, yr=yr, n_pix=n_pix, mode=mode, kernel=kernel, bandwidth=bandwidth, bintype=bintype
	;;-----
	;; Get 2d map from particles
	;;
	;;	mode: [1] integer
	;;		positive (>0) gives denstiy at each point
	;;		negative (<0) gives density at each grid
	;;
	;;	kernel: [1] integer
	;;		0 - 2D histogram
	;;		1 - Gaussian kernel
	;;		2 - Triangular kernel (not implemented)
	;;		3 - p^3M (not implemented)
	;;
	;;	bintype: [1] string
	;;		binning scheme
	;;		'NGP' or 'CIC'
	;;
	;;	bandwidth: [2] single / double
	;;		kernel size
	;;		if not set, a pixel size used
	;;-----



	;;-----
	;; Initial settings
	;;-----
	IF ~KEYWORD_SET(zz) THEN zz = xx*0.d + 1.d
	IF ~KEYWORD_SET(n_pix) THEN n_pix = 1000L
	IF ~KEYWORD_SET(mode) THEN mode = -1L
	IF ~KEYWORD_SET(kernel) THEN kernel = 1L
	IF ~KEYWORD_SET(bandwidth) THEN bandwidth = [xr(1)-xr(0), yr(1)-yr(0)]/n_pix
	IF ~KEYWORD_SET(bintype) THEN bintype = 'NGP'

	xx = DOUBLE(xx) & yy = DOUBLE(yy) & zz = DOUBLE(zz)
	xr = DOUBLE(xr) & yr = DOUBLE(yr) & n_pix = LONG(n_pix)

	n_dim = 2L
	hmat 	= DBLARR(2,2)
	hmat(0,0)	= bandwidth(0)^2
	hmat(1,1) 	= bandwidth(1)^2
	ihmat 	= INVERT(hmat)
	const 	= DOUBLE(DETERM(hmat)^(-0.5))*(!pi*2.)^(-n_dim/2.)

	settings	= self->getheader()

	;;-----
	;; Range cut
	;;-----
	cut 	= WHERE(xx GE xr(0) AND xx LE xr(1) AND yy GE yr(0) AND yy LE yr(1), ncut)
	IF ncut EQ 0L THEN BEGIN
		self->errorout, 'NO DATA ARE LEFT WITH THE GIVEN RANGES'
		RETURN, -1
	ENDIF

	dx = xx(cut) & dy = yy(cut) & dz = zz(cut)

	;;-----
	;; Compute
	;;-----
	density 	= FLTARR(n_pix, n_pix)
	;ptcl 		= FLTARR(n_pix)

	IF KERNEL EQ 0L THEN BEGIN ;; 2D histogram
		ftr_name 	= settings.dir_lib + '/src/fortran/js_kde_2dhisto.so'
			larr = LONARR(20)
			larr(0) = N_ELEMENTS(dx)
			larr(1) = settings.num_thread
			larr(2) = n_pix

		void = CALL_EXTERNAL(ftr_name, 'js_kde_2dhisto', dx, dy, dz, density, ptcl, xr, yr, larr)
	ENDIF ELSE IF kernel EQ 1L THEN BEGIN ;; Gaussian Kernel

		;;----- Binning
		grid2 	= DBLARR(n_pix*2, n_pix*2)

		ftr_name 	= settings.dir_lib + '/src/fortran/js_kde_gauss_binning.so'
			larr= LONARR(20)
			larr(0) = N_ELEMENTS(dx)
			larr(1) = settings.num_thread
			larr(2) = n_pix

			IF bintype EQ 'NGP' THEN larr(3) = 1L
			IF bintype EQ 'CIC' THEN larr(3) = 2L

			larr(4) = 1L

		void 	= CALL_EXTERNAL(ftr_name, 'js_kde_gauss_binning', $
			dx, dy, dz, xr, yr, grid2, larr)

		gridorg 	= grid2

		;;----- Kernel Matrix
		delX 	= (xr(1)-xr(0))/n_pix
		delY 	= (yr(1)-yr(0))/n_pix

		grid_ker	= DBLARR(n_pix*2., n_pix*2.)

		ix 	= DINDGEN(n_pix*2) + 0.5
		iy 	= DINDGEN(n_pix*2) + 0.5

		ix 	= REBIN(ix, n_pix*2, n_pix*2)
		iy 	= REBIN(TRANSPOSE(iy), n_pix*2, n_pix*2)

		ix 	= (n_pix - 0.5)*delX - ix * delX
		iy 	= (n_pix - 0.5)*delY - iy * delY

		dum_ker 	= -0.5*(ihmat(0,0) * ix^2 + ihmat(1,1) * iy^2 + $
			(ihmat(0,1) + ihmat(1,0))*ABS(ix)*ABS(iy))

		avoid 	= WHERE(dum_ker GT -700.0, nn)
		IF nn EQ 0L THEN BEGIN 
			PRINT, '?' 
			STOP
		ENDIF

		grid_ker(avoid) 	= EXP(dum_ker(avoid))*const

		;;----- Convolution
		ift_bin 	= FFT(grid2, /DOUBLE)
		ift_ker		= FFT(grid_ker, /DOUBLE)

		ift_conv 	= ift_bin * ift_ker
		grid 	= FFT(ift_conv, /DOUBLE, /inverse)
		grid 	= DOUBLE(REAL_PART(grid))
		grid2 = 0. & grid_ker = 0. & ift_bin = 0. & ift_ker = 0. & ix = 0. & iy = 0.

		;;----- Adjust
		grid2 	= grid
		grid(*, 0:n_pix-1L) 	= grid2(*,n_pix:*)
		grid(*,n_pix:*)			= grid2(*,0:n_pix-1L)

		grid2 	= grid
		grid(0:n_pix-1L,*) 		= grid2(n_pix:*,*)
		grid(n_pix:*,*)			= grid2(0:n_pix-1L,*)

		grid 	= grid(n_pix/2-1:n_pix/2+n_pix-2,n_pix/2-1:n_pix/2+n_pix-2)
		density 	= grid

		grid 	= 0.

		IF mode GT 0. THEN BEGIN
			PRINT, 'interpolation'

			ptcl 	= DBLARR(N_ELEMENTS(dx))

			ftr_name	= settings.dir_lib + '/src/fortran/js_kde_gauss_pts.so'
				larr= LONARR(20)
				larr(0) = N_ELEMENTS(dx)
				larr(1) = settings.num_thread
				larr(2) = n_pix

			void 	= CALL_EXTERNAL(ftr_name, 'js_kde_gauss_pts', $
				dx, dy, dz, xr, yr, density, ptcl, larr)
		ENDIF
	ENDIF

	;;-----
	;; RETURN
	;;-----

	IF mode GT 0L THEN BEGIN
		ptcl 	= ptcl / TOTAL(ptcl) * TOTAL(dz)
		RETURN, {x:dx, y:dy, z:ptcl}
	ENDIF ELSE BEGIN
		delx 	= (xr(1)-xr(0))/n_pix
		dely 	= (yr(1)-yr(0))/n_pix

		ix 	= FINDGEN(n_pix) + 0.5
		iy 	= FINDGEN(n_pix) + 0.5
		ix 	= REBIN(ix, n_pix, n_pix)
		iy 	= REBIN(TRANSPOSE(iy), n_pix, n_pix)
		ix 	= ix*delx + xr(0)
		iy 	= iy*dely + yr(0)

		cut_neg 	= WHERE(density LT 0., ncut)
		IF ncut GE 1L THEN density(cut_neg) = 0.d

		density 	= density / TOTAL(density) * TOTAL(dz) / delx / dely

		RETURN, {x:ix, y:iy, z:density}
	ENDELSE
END

FUNCTION veluga::d_gasmap, n_snap, cell, xr, yr, n_pix=n_pix, $
	amrvar=amrvar, amrtype=amrtype, minlev=minlev, maxlev=maxlev, proj=proj, $
	delZ=delZ, xx0=xx0, vv0=vv0, memeff=memeff
	;celltype=celltype, cellphase=cellphase, 

	;;-----
	;; Get gas map from amr data
	;;	cell data should be given in physical unit
	;;
	;;	amrvar: [1] string
	;;		'D' 	- density
	;;		'T'		- temperature
	;;		'PT'	- Thermal presusre
	;;		'PR'	- Extenral pressure by rho X v^2 (xx0 and vv0 should be argued)
	;;		'Z'		- Metallicity
	;;
	;;	amrtype: [1] string
	;;		'MW'	- mass weighted
	;;		'VW'	- volume weighted
	;;		'MAX'	- Maximum along LOS
	;;		'CD'	- Column Denstiy
	;;		'HIST'	- Histogram
	;;
	;;	minlev, maxlev: [1] long
	;;		min and max amr level (if not set, a simulation value is employed)
	;;
	;;	proj: [1] string
	;;		'xy'	- xy
	;;
	;;	delZ: [1] double
	;;		slicing width
	;;
	;;	xx0, vv0: [3] double
	;;		central position and velocity for computing external pressure
	;;
	;;	memeff: boolean
	;;		if argued, export the original quantitiy map and density map
	;;-----

	;;-----
	;; Initial settings
	;;-----

	info  	= self->g_info(n_snap)
	settings= self->getheader()

	IF ~KEYWORD_SET(n_pix) THEN n_pix = 1000L
	IF ~KEYWORD_SET(amrvar) THEN amrvar = 'D'
	IF ~KEYWORD_SET(amrtype) THEN amrtype = 'MW'
	IF ~KEYWORD_SET(delZ) THEN delZ = -1.d

	

	IF ~KEYWORD_SET(minlev) THEN minlev = info.levmin
	IF ~KEYWORD_SET(maxlev) THEN maxlev = info.levmax

	amrvar	= [STRUPCASE(amrvar)]
	amrtype = [STRUPCASE(amrtype)]
	amrtype_l	= LONARR(N_ELEMENTS(amrvar))

	FOR ai=0L, N_ELEMENTS(amrvar)-1L DO BEGIN
		CASE amrtype(ai) OF
			'MW'  : amrtype_l(ai) = 1
			'VW'  : amrtype_l(ai) = 2
			'MAX' : amrtype_l(ai) = 3
			'CD'  : amrtype_l(ai) = 4
			'HIST': amrtype_l(ai) = 5
		ENDCASE
	ENDFOR

	cut 	= WHERE(amrtype EQ 'PR', ncut)
	IF ncut GE 1L THEN BEGIN
		IF ~KEYWORD_SET(xx0) THEN BEGIN
			xx0 = [0.d, 0.d, 0.d]
			self->errourout, 'Centarl position is not given to compute extenral pressure'
		ENDIF
		IF ~KEYWORD_SET(vv0) THEN BEGIN
			vv0 = [0.d, 0.d, 0.d]
			self->errourout, 'Centarl velocity is not given to compute extenral pressure'
		ENDIF
	ENDIF

	IF delz GT 0. AND ~KEYWORD_SET(xx0) THEN BEGIN
		self->errorout, 'Cental position should be given to determine thickness along LOS'
		xx0 	= [0.d, 0.d, 0.d]
	ENDIF


	;;-----
	;; ALLOCATE
	;;-----
	temp 	= DBLARR(cell.n, N_ELEMENTS(amrvar), 2)

	;; nCell X nAMR X [weight, map] 
	;temp 	= DBLARR(cell.N, 2)
		;; 0 as variable
		;; 1 as density (for MW)

	;temp(*,0) 	= cell.den

	;; Sun chemistry (Asplund 09)
	sun_N_o_H	= 10.d^(7.83d - 12.d) ;; [#_N / #_H]
	sun_C_o_H	= 10.d^(8.43d - 12.d)
	sun_O_o_H	= 10.d^(8.69d - 12.d)
	sun_Mg_o_H	= 10.d^(7.60d - 12.d)
	sun_Si_o_H	= 10.d^(7.51d - 12.d)
	sun_S_o_H	= 10.d^(7.12d - 12.d)
	sun_Fe_o_H	= 10.d^(7.50d - 12.d)

	FOR i=0L, N_ELEMENTS(amrvar)-1L DO BEGIN
		CASE amrvar(i) OF
			'D' : temp(*,i,0) = cell.den
			'T' : temp(*,i,0) = cell.temp
			'PT': temp(*,i,0) = cell.P_thermal
			'PR': temp(*,i,0) = cell.den * (self->g_d3d(cell.vx, cell.vy, cell.vz, vv0))^2
			'Z' : temp(*,i,0) = cell.zp
			'O/FE' 	: BEGIN
				N_X 	= cell.chem_O / 15.999d
				N_H 	= cell.chem_H
				N_Fe 	= cell.chem_Fe / 55.845d
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(N_Fe / N_H) - (ALOG10(sun_O_o_H) - ALOG10(sun_Fe_o_H))
				END
			'MG/FE' 	: BEGIN
				N_X 	= cell.chem_Mg / 24.305d
				N_H 	= cell.chem_H
				N_Fe 	= cell.chem_Fe / 55.845d
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(N_Fe / N_H) - (ALOG10(sun_Mg_o_H) - ALOG10(sun_Fe_o_H))
				END
			'SI/FE' 	: BEGIN
				N_X 	= cell.chem_Si / 28.0855d
				N_H 	= cell.chem_H
				N_Fe 	= cell.chem_Fe / 55.845d
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(N_Fe / N_H) - (ALOG10(sun_Si_o_H) - ALOG10(sun_Fe_o_H))
				END
			'ALPHA/FE'	: BEGIN
				N_X 	= cell.chem_Si / 28.0855d + cell.chem_O / 15.999d + cell.chem_Mg / 24.305d
				N_H 	= cell.chem_H
				N_Fe 	= cell.chem_Fe / 55.845d * 3.d
				sOh 	= (sun_Mg_o_H + sun_Si_o_H + sun_O_o_H)
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(N_Fe / N_H) - (ALOG10(sOh) - ALOG10(sun_Fe_o_H*3.d))
				END
			'S/FE' 	: BEGIN
				N_X 	= cell.chem_S / 32.065d
				N_H 	= cell.chem_H
				N_Fe 	= cell.chem_Fe / 55.845d
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(N_Fe / N_H) - (ALOG10(sun_S_o_H) - ALOG10(sun_Fe_o_H))
				END
			'C/H' 	: BEGIN
				N_X 	= cell.chem_C / 12.011d
				N_H 	= cell.chem_H
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(sun_C_o_H)
				END
			'N/H' 	: BEGIN
				N_X 	= cell.chem_N / 14.0067d
				N_H 	= cell.chem_H
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(sun_N_o_H)
				END
			'LIGHT/H' 	: BEGIN
				N_X 	= cell.chem_N / 14.0067d + cell.chem_C / 12.011d
				N_H 	= cell.chem_H * 2.d
				soH 	= sun_N_o_H + sun_C_o_H
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(soh/2.)
				END
			'FE/H' 	: BEGIN
				N_X 	= cell.chem_Fe / 55.845d
				N_H 	= cell.chem_H
				temp(*,i,0) 	= ALOG10(N_X / N_H) - ALOG10(sun_Fe_o_H)
				END
			'O'		: temp(*,i,0)	= cell.chem_O * cell.mp / 15.999d
			'SI'	: temp(*,i,0)	= cell.chem_Si * cell.mp / 28.0855d
			'MG'	: temp(*,i,0)	= cell.chem_Mg * cell.mp / 24.305d
			'FE'	: temp(*,i,0)	= cell.chem_Fe * cell.mp / 55.845d
	
			'DUST1' : temp(*,i,0)	= cell.mp * cell.dust_1
			'DUST2' : temp(*,i,0)	= cell.mp * cell.dust_2
			'DUST3' : temp(*,i,0)	= cell.mp * cell.dust_3
			'DUST4' : temp(*,i,0)	= cell.mp * cell.dust_4
			'DUST' 	: temp(*,i,0)	= cell.mp * (cell.dust_1 + cell.dust_2 + cell.dust_3 + cell.dust_4)
		ENDCASE

		CASE amrtype(i) OF
			'MW'	: temp(*,i,1) 	= cell.den 	;; converted to mass
			ELSE 	: temp(*,i,1) 	= cell.den 		;; actually not used
		ENDCASE
	ENDFOR

	map 	= DBLARR(n_pix, n_pix, N_ELEMENTS(amrvar), 2)

	;;-----
	;; Compute
	;;-----
	levind 	= cell.levelind

	integrity 	= 0L
	FOR lev=minlev, maxlev DO BEGIN
		IF levind(lev,2) EQ 0L THEN CONTINUE
		ind0 	= levind(lev,0)
		ind1 	= levind(lev,1)

		xx 	= cell.xx(ind0:ind1)
		yy 	= cell.yy(ind0:ind1)
		zz 	= cell.zz(ind0:ind1)
		level= cell.level(ind0:ind1)
		tempdum		= temp(ind0:ind1,*,*)
		dx 	= cell.dx(ind0)

		check 	= ABS(level - lev)
		IF MAX(check) GT 0L THEN BEGIN
			cut 	= WHERE(cell.level EQ lev, nlev)
			IF nlev EQ 0L THEN CONTINUE
			xx 	= cell.xx(cut)
			yy 	= cell.yy(cut)
			zz 	= cell.zz(cut)
			tempdum	= temp(cut,*,*)
			dx 	= cell.dx(cut(0))
			integrity	+= nlev
			print, 'here?'
		ENDIF ELSE BEGIN
			integrity	+= levind(lev,2)
		ENDELSE

		IF delz GT 0. THEN BEGIN
			CASE proj OF
				'xy' : dz = ABS(zz - xx0(2))
				'xz' : dz = ABS(yy - xx0(1))
				'yz' : dz = ABS(xx - xx0(0))
			ENDCASE

			cut 	= WHERE(dz LT delz, nlev)
			
			xx 	= xx(cut)
			yy 	= yy(cut)
			zz 	= zz(cut)
			tempdum	= tempdum(cut,*,*)
		ENDIF

		;IF nlev EQ 0L THEN CONTINUE

		bandwidth 	= [1.d, 1.d]*dx

		ftr_name 	= settings.dir_lib + '/src/fortran/js_gasmap.so'
		larr = LONARR(20) & darr = DBLARR(20)

		larr(0)	= N_ELEMENTS(xx)
		larr(1) = N_ELEMENTS(amrvar)
		larr(2) = n_pix
		larr(3) = self.num_thread

		darr(0) = info.cgs.kpc / ((xr(1)-xr(0))/n_pix*(yr(1)-yr(0))/n_pix)
		;; column density unit conversion (> /cm^2)

		CASE proj OF
			'xy' : BEGIN
				xx2 = xx
				yy2 = yy
				END
			'xz' : BEGIN
				xx2 = xx
				yy2 = zz
				END
			'yz' : BEGIN
				xx2 = yy
				yy2 = zz
				END
			ELSE: BEGIN
				self->errourout,'Proper proj should be given: stop here'
				STOP
				END
		ENDCASE

		

		void 	= CALL_EXTERNAL(ftr_name, 'js_gasmap', $
			larr, darr, xx2, yy2, tempdum, bandwidth, DOUBLE(xr), DOUBLE(yr), map, amrtype_l)

	ENDFOR

	IF integrity NE N_ELEMENTS(temp(*,0,0)) THEN BEGIN
		self->errorout, 'levelind integrity is broken'
		self->errorout, 'N_cell = ', STRTRIM(N_ELEMENTS(temp(*,0,0)),2)
		self->errorout, 'N_lev  = ', STRTRIM(integrity,2)
		STOP
	ENDIF

	;denmap 	= REFORM(map(*,*,0), n_pix, n_pix)
	;map 	= REFORM(map(*,*,0), n_pix, n_pix)

	;;----- output
	;dummymap	= DBLARR(n_pix, n_pix)
	result 	= REPLICATE({amrvar:'', amrtype:'', map:DBLARR(n_pix,n_pix), map0:DBLARR(n_pix,n_pix)}, N_ELEMENTS(amrvar))

	FOR i=0L, N_ELEMENTS(amrvar)-1L DO BEGIN
		result(i).amrvar = amrvar(i)
		result(i).amrtype= amrtype(i)
		result(i).map 	= REFORM(map(*,*,i,0),n_pix, n_pix)
		result(i).map0 	= REFORM(map(*,*,i,1),n_pix, n_pix)
	ENDFOR
	
	IF KEYWORD_SET(memeff) THEN RETURN, result

	FOR i=0L, N_ELEMENTS(amrtype)-1L DO BEGIN
		cut 	= WHERE(result(i).map0 GT 0., ncut)
		IF ncut EQ 0L THEN CONTINUE

		CASE amrtype(i) OF
			'MW': result(i).map(cut) /= result(i).map0(cut)
			'VW': result(i).map(cut) /= result(i).map0(cut)
			'MAX':
			'CD':
			'HIST':
		ENDCASE
	ENDFOR

	FOR i=0L, N_ELEMENTS(amrtype)-1L DO BEGIN
		cut 	= WHERE(result(i).map0 EQ 0., ncut)
		IF ncut EQ 0L THEN CONTINUE

		CASE amrtype(i) OF
			'MW': result(i).map(cut) = 0.d
			'VW': result(i).map(cut) = 0.d
			'MAX':
			'CD':
			'HIST':
		ENDCASE
	ENDFOR

	RETURN, result
END

;;-----
;; TABLE GENERATOR & LOAD
;;-----
FUNCTION t_conformal_sub, A, _extra=extra
	oM 	= extra.oM
	oL 	= extra.oL
	RETURN, 1./(A^3. * SQRT(oM / A^3. + oL))
END
PRO veluga::t_conformal, oM, oL
	settings 	= self->getheader()
	sfact	= DINDGEN(10000)/9999.*0.98 + 0.02 
	conft	= DBLARR(10000)

	FOR i=0L, N_ELEMENTS(sfact)-1L DO BEGIN
		QSIMP, 't_conformal_sub', sfact(i), 1, val, /double, oM=oM, oL=oL
		conft(i) 	= val * (-1.)
	ENDFOR

	soM 	= STRING(oM*1000., format='(I3.3)')
	soL 	= STRING(oL*1000., format='(I3.3)')
	fname  	= settings.dir_lib + 'table/conformal_table_' + soM + '_' + soL + '.sav'
	
	tbl 	= {sfact:sfact, conft:conft, oM:oM, oL:oL}
	SAVE, filename=fname, tbl
END
FUNCTION veluga::t_conformal_load, oM, oL
	settings 	= self->getheader()
	soM 	= STRING(oM*1000., format='(I3.3)')
	soL 	= STRING(oL*1000., format='(I3.3)')
	fname  	= settings.dir_lib + 'table/conformal_table_' + soM + '_' + soL + '.sav'

	isfile 	= FILE_SEARCH(fname)
	IF STRLEN(isfile) LE 5L THEN BEGIN
		self->t_conformal, oM, oL
	ENDIF

	RESTORE, fname
	RETURN, tbl
END


;;

FUNCTION t_lbt_sub, X, _extra=extra
	oM	= extra.oM
	oL	= extra.oL

	RETURN, 1./(1.+X)/SQRT(OM*(1.+X)^3 + OL)
END
PRO veluga::t_lbt, oM, oL, H0
	settings 	= self->getheader()

	tmp_red = DINDGEN(10000)/9999.*0.98 + 0.02
	tmp_red = 1./tmp_red - 1. & tmp_red = REVERSE(tmp_red)
	tmp_gyr = tmp_red & tmp_gyr(0) = 0.

	FOR i=1L, N_ELEMENTS(tmp_red)-1L DO BEGIN
		QSIMP, 't_lbt_sub', 0., tmp_red(i), val, oM=oM, oL=oL
		tmp_gyr(i) 	= val
	ENDFOR
	tmp_red(0) 	= 0.d
	tmp_gyr(0) 	= 0.d
	tmp_gyr	= tmp_gyr / H0 *3.08568025e19 / 3.1536000d+16

	tbl 	= {redsh:tmp_red, gyr:tmp_gyr, oM:oM, oL:oL, H0:H0}
	
	soM 	= STRING(oM*1000., format='(I3.3)')
	soL 	= STRING(oL*1000., format='(I3.3)')
	sH0		= STRING(H0*10., format='(I3.3)')
	fname  	= settings.dir_lib + 'table/lbt_table_' + soM + '_' + soL + '_' + sH0 + '.sav'
	SAVE, filename=fname, tbl
END
FUNCTION veluga::t_lbt_load, oM, oL, H0
	settings	= self->getheader()
	soM 	= STRING(oM*1000., format='(I3.3)')
	soL 	= STRING(oL*1000., format='(I3.3)')
	sH0		= STRING(H0*10., format='(I3.3)')

	fname  	= settings.dir_lib + 'table/lbt_table_' + soM + '_' + soL + '_' + sH0 + '.sav'

	isfile 	= FILE_SEARCH(fname)
	IF STRLEN(isfile) LE 5L THEN BEGIN
		self->t_lbt, oM, oL, H0
	ENDIF

	RESTORE, fname
	RETURN, tbl
END

;;
FUNCTION t_comv_sub, X, _extra=extra
	oM 	= extra.oM
	oL 	= extra.oL

	RETURN, 1./SQRT(oM * (1+X)^3 + oL)
END
PRO veluga::t_comv, oM, oL, H0
	settings 	= self->getheader()

	tmp_red = DINDGEN(10000)/9999.*0.98 + 0.02
	tmp_red = 1./tmp_red - 1. & tmp_red = REVERSE(tmp_red)
	tmp_comv 	= tmp_red

	FOR i=1L, N_ELEMENTS(tmp_red)-1L DO BEGIN
		QSIMP, 't_comv_sub', 0., tmp_red(i), val, oM=oM, oL=oL
		tmp_comv(i) 	= val
	ENDFOR
	tmp_comv(0)	= 0.d
	tmp_red(0)	= 0.d
	tmp_comv 	= tmp_comv * 2.99792458e5 / H0

	tbl 	= {redsh:tmp_red, comv:tmp_comv, oM:oM, oL:oL, H0:H0}
	soM 	= STRING(oM*1000., format='(I3.3)')
	soL 	= STRING(oL*1000., format='(I3.3)')
	sH0		= STRING(H0*10., format='(I3.3)')
	fname  	= settings.dir_lib + 'table/comv_table_' + soM + '_' + soL + '_' + sH0 + '.sav'
	SAVE, filename=fname, tbl
END
FUNCTION veluga::t_comv_load, oM, oL, H0
	settings	= self->getheader()
	soM 	= STRING(oM*1000., format='(I3.3)')
	soL 	= STRING(oL*1000., format='(I3.3)')
	sH0		= STRING(H0*10., format='(I3.3)')

	fname  	= settings.dir_lib + 'table/comv_table_' + soM + '_' + soL + '_' + sH0 + '.sav'

	isfile 	= FILE_SEARCH(fname)
	IF STRLEN(isfile) LE 5L THEN BEGIN
		self->t_comv, oM, oL, H0
	ENDIF

	RESTORE, fname
	RETURN, tbl
END

;;

PRO veluga::t_miles_sdss
	settings 	= self->getheader()

	fname 	= settings.dir_lib + 'table/sdss_ch_iPp0.00.MAG'
	READCOL, fname, dum1, dum2, metal, age, mass, u, g, r, i, z, u2, g2, r2, i2, z2, format='A, D, D, D, D, D, D, D, D, D, D, D, D, D, D', skipline=1L

	tbl 	= REPLICATE({age:0.d, metal:0.d, mass:0.d, u:0.d, g:0.d, r:0.d, i:0.d, z:0.d}, N_ELEMENTS(metal))
	tbl.age 	= age
	tbl.metal 	= (10.^metal)*settings.sun_met
	tbl.u 		= u2
	tbl.g 		= g2
	tbl.r 		= r2
	tbl.i 		= i2
	tbl.z 		= z2

	age_arr 	= self->g_unique(tbl.age)
	met_arr 	= self->g_unique(tbl.metal)

	;; Make 2D TABLE

	ind 	= SORT(tbl.metal)
	tbl 	= tbl(ind)

	n_met 	= N_ELEMENTS(met_arr)
	n_age 	= N_ELEMENTS(age_arr)

	u 	= REFORM(tbl.u, n_age, n_met)
	g 	= REFORM(tbl.g, n_age, n_met)
	r 	= REFORM(tbl.r, n_age, n_met)
	i 	= REFORM(tbl.i, n_age, n_met)
	z 	= REFORM(tbl.z, n_age, n_met)

	FOR k=0L, n_met-1L DO BEGIN
		dum 	= u(*,k)
		cut 	= WHERE(dum LT 0., ncut)
		IF ncut GE 1L THEN BEGIN
			dum(cut)	= MIN(dum(WHERE(dum GT 0.)))
			u(*,k) 	= dum
		ENDIF

		dum 	= g(*,k)
		cut 	= WHERE(dum LT 0., ncut)
		IF ncut GE 1L THEN BEGIN
			dum(cut)	= MIN(dum(WHERE(dum GT 0.)))
			g(*,k) 	= dum
		ENDIF

		dum 	= r(*,k)
		cut 	= WHERE(dum LT 0., ncut)
		IF ncut GE 1L THEN BEGIN
			dum(cut)	= MIN(dum(WHERE(dum GT 0.)))
			r(*,k) 	= dum
		ENDIF

		dum 	= i(*,k)
		cut 	= WHERE(dum LT 0., ncut)
		IF ncut GE 1L THEN BEGIN
			dum(cut)	= MIN(dum(WHERE(dum GT 0.)))
			i(*,k) 	= dum
		ENDIF

		dum 	= z(*,k)
		cut 	= WHERE(dum LT 0., ncut)
		IF ncut GE 1L THEN BEGIN
			dum(cut)	= MIN(dum(WHERE(dum GT 0.)))
			z(*,k) 	= dum
		ENDIF
	ENDFOR

	ref 	= {age:age_arr, metal:met_arr, u:u, g:g, r:r, i:i, z:z}
	SAVE, filename=settings.dir_lib + 'table/miles_sdss.sav', ref
	RETURN
END
FUNCTION veluga::t_miles_sdss_load
	settings = self->getheader()
	fname 	= settings.dir_lib + 'table/miles_sdss.sav'
	isfile 	= FILE_SEARCH(fname)
	IF STRLEN(isfile) LE 5L THEN BEGIN
		self->t_miles_sdss
	ENDIF

	RESTORE, fname
	RETURN, ref
END

PRO veluga::t_miles_galex
	settings 	= self->getheader()

	fname 	= settings.dir_lib + 'table/NUV_ch_iPp0.00.MAG'
	READCOL, fname, dum1, dum2, metal, age, mass, nuv, nuv2, format='A, D, D, D, D, D, D', skipline=1L

	tbl 	= REPLICATE({age:0.d, metal:0.d, mass:0.d, nuv:0.d}, N_ELEMENTS(metal))

	tbl.age 	= age
	tbl.metal 	= (10.^metal)*settings.sun_met
	tbl.nuv 	= nuv2

	age_arr		= self->g_unique(tbl.age)
	met_arr		= self->g_unique(tbl.metal)

	;; Make 2D TABLE

	ind 	= SORT(tbl.metal)
	tbl 	= tbl(ind)

	n_met 	= N_ELEMENTS(met_arr)
	n_age 	= N_ELEMENTS(age_arr)

	nuv 	= REFORM(tbl.nuv, n_age, n_met)
	FOR i=0L, n_met-1L DO BEGIN
		dum 	= nuv(*,i)
		cut 	= WHERE(dum LT 0., ncut)
		IF ncut GE 1L THEN BEGIN
			dum(cut)	= MIN(dum(WHERE(dum GT 0.)))
			nuv(*,i) 	= dum
		ENDIF
	ENDFOR

	ref 	= {age:age_arr, metal:met_arr, nuv:nuv}
	SAVE, filename=settings.dir_lib + 'table/miles_galex.sav', ref
	RETURN
END
FUNCTION veluga::t_miles_galex_load
	settings = self->getheader()
	fname 	= settings.dir_lib + 'table/miles_galex.sav'
	isfile 	= FILE_SEARCH(fname)
	IF STRLEN(isfile) LE 5L THEN BEGIN
		self->t_miles_galex
	ENDIF

	RESTORE, fname
	RETURN, ref
END
	

;;-----
;; MAIN
;;-----
PRO veluga__define

	void	= {veluga, header:PTR_NEW(), fname:PTR_NEW(), tree:PTR_NEW(), num_thread:1L}

	;self->rdheader, *self,fname
	;veluga->rdheader, *veluga.fname

	RETURN
END

;VELociraptor Utilities for Galaxy Analysis
