FUNCTION veluga::init, fname, num_thread=num_thread

	IF ~KEYWORD_SET(num_thread) THEN num_thread = 1L

	self.header	= PTR_NEW(/allocate)
	self.rdheader, fname
	self.num_thread	= num_thread

	ptr 	= REPLICATE({tree:PTR_NEW(), key:PTR_NEW(), stat:-2L},2)
	self.tree 	= PTR_NEW(ptr)
	
	settings 	= self.getheader()

	RETURN, 1
END

FUNCTION veluga::allocate, nn, type=type

	CASE type OF
		'part'		: RETURN, REPLICATE({xx:0.d, yy:0.d, zz:0.d, vx:0.d, vy:0.d, vz:0.d, mp:0.d, ap:0.d, zp:0.d, id:0L, family:0L, domain:0L, KE:0.d, UE:0.d, PE:0.d}, nn)
		'cell'		: RETURN, REPLICATE({xx:0.d, yy:0.d, zz:0.d, vx:0.d, vy:0.d, vz:0.d, level:0L, dx:0.d, den:0.d, temp:0.d, zp:0.d, mp:0.d, KE:0.d, UE:0.d, PE:0.d}, nn)
		ELSE: STOP
	ENDCASE
END
;;-----
;; HEADER
;;-----
PRO veluga::rdheader, fname

	FINDPRO, 'veluga__define', dirlist=dirlist, /noprint
	
	settings 	= {header:'^^', dir_lib:dirlist(0), $
		erg_msg_i:'%123123----- VELUGA -----', erg_msg_f:'%123123----------', erg_msg_0:' '}

	
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
;;-----
;; READ ROUTINE
;;-----
FUNCTION veluga::r_gal, snap0, id0, horg=horg

	;;-----
	;; READ Galaxies
	;;-----
	IF ~KEYWORD_SET(horg) THEN horg = 'g'
	settings	= self.getheader()
	dir 		= settings.dir_catalog

	;;-----
	;; LOAD COLUMN
	;;-----
	gprop		= [settings.column_list, settings.gal_prop]

	FOR i=0L, N_ELEMENTS(Gprop)-1L DO BEGIN
		IF Gprop(i) EQ 'sfr' THEN Gprop(i) = 'SFR'
		IF Gprop(i) EQ 'ABMAG' OR $
			Gprop(i) EQ 'abmag' THEN Gprop(i) = 'ABmag'
	ENDFOR

	IF horg EQ 'g' THEN $
		fname = dir + 'Galaxy/VR_Galaxy/snap_' + STRING(snap0,format='(I4.4)') + '.hdf5'
	IF horg EQ 'h' THEN $
		fname = dir + 'Halo/VR_Halo/snap_' + STRING(snap0,format='(I4.4)') + '.hdf5'

	;;-----
	;; OPEN HDF5
	;;-----
	fid	= H5F_OPEN(fname)

	;;-----
	;; ID AND MASS
	;;-----
	did	= H5D_OPEN(fid, 'Mass_tot')
	mass_tot 	= H5D_READ(did)
	H5D_CLOSE, did

	did	= H5D_OPEN(fid, 'Mvir')
	mvir 	= H5D_READ(did)
	H5D_CLOSE, did

	did	= H5D_OPEN(fid, 'ID')
	ID 	= H5D_READ(did)
	H5D_CLOSE, did

	IF id0 LT 0L THEN BEGIN 
		IF KEYWORD_SET(masscut) THEN BEGIN
			IF horg EQ 'g' THEN $
				mcut	= WHERE(mass_tot GE masscut(0) AND mass_tot LT masscut(1), nn)
			IF horg EQ 'h' THEN $
				mcut	= WHERE(mvir GE masscut(0) AND mvir LT masscut(1), nn)

			IF nn EQ 0L THEN BEGIN
				PRINT, '%123123-----'
				PRINT, '	NO GALAXIES (HALOS) WITHIN THE MASS LIMIT'
				PRINT, '	(CONVERT TO READING ALL GALS'

				mcut	= WHERE(mvir GE 0. AND mvir LT 1e20, nn)
			ENDIF
			n_gal	= nn
		ENDIF ELSE BEGIN
			mcut	= WHERE(mvir GE 0. AND mvir LT 1e20, nn)
			n_gal	= nn
		ENDELSE
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

	did	= H5D_OPEN(fid, 'ID_000001/Domain_List')
	dlist	= H5D_READ(did)
	H5D_CLOSE, did

	did	= H5D_OPEN(fid, 'CONF_R')
	conf_r	= H5D_READ(did)
	H5D_CLOSE, did

	tmpstr	+= 'rate:1.0d, Domain_List:dlist, Aexp:1.0d, CONF_R:CONF_R'
	n_mpi	= N_ELEMENTS(dlist)

	did	= H5D_OPEN(fid, 'Flux_List')
	flux_list	= H5D_READ(did)
	H5D_CLOSE, did

	did	= H5D_OPEN(fid, 'MAG_R')
	mag_r	= H5D_READ(did)
	H5D_CLOSE, did

	did	= H5D_OPEN(fid, 'SFR_R')
	sfr_r	= H5D_READ(did)
	H5D_CLOSE, did

	did	= H5D_OPEN(fid, 'SFR_T')
	sfr_t	= H5D_READ(did)
	H5D_CLOSE, did

	tmpstr	+= ', isclump:-1L, Flux_List:flux_list, '
	tmpstr	+= 'MAG_R:MAG_R, SFR_R:SFR_R, SFR_T:SFR_T}'

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

		did	= H5D_OPEN(fid, idstr + '/isclump')
		GP(i).isclump	= H5D_READ(did)
		H5D_CLOSE, did


		did	= H5D_open(fid, 'Flux_List')
		GP(i).flux_list	= H5D_READ(did)
		H5D_close, did

		did	= H5D_OPEN(fid, idstr + '/Domain_List')
		GP(i).Domain_list	= H5D_READ(did)
		H5D_CLOSE, did

		did	= H5D_OPEN(fid, idstr + '/rate')
		GP(i).rate	= H5D_READ(did)
		H5D_CLOSE, did

		did	= H5D_open(fid, idstr + '/Aexp')
		GP(i).aexp	= H5D_READ(did)
		H5D_close, did

		GP(i).snapnum	= snap0
		GP(i).redsh 	= 1./GP(i).aexp - 1.d
	ENDFOR
	H5F_CLOSE, fid
	RETURN, GP
END

FUNCTION veluga::r_ptcl, snap0, id0, horg=horg, simout=simout

	;;-----
	;; READ MEMBER PTCL
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

	ftr_name 	= settings.dir_lib + 'fortran/get_ptcl.so'
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
		self->errorout, '		r_ptcl: NO MATCHED PTCLs'
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
		output.xx 	*= (info.unit_l/3.086e21)
		output.yy 	*= (info.unit_l/3.086e21)
		output.zz 	*= (info.unit_l/3.086e21)

		output.vx 	*= (info.kms)
		output.vy 	*= (info.kms)
		output.vz 	*= (info.kms)

		output.mp 	*= (info.unit_m / 1.98892d33)
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
		self->errout, 'No tree data exists for ' + horg
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
;;-----
;; SIMPLE GET FTNS
;;	-- RAMSES RELATED
;;-----
FUNCTION veluga::g_info, snap0
	str 	= STRING(snap0, format='(I5.5)')
	settings 	= self.getheader()
	rd_info, info, file=settings.dir_raw + 'output_' + str + '/info_' + str + '.txt'
	RETURN, info
END

FUNCTION veluga::g_domain, xc2, yc2, zc2, rr2, snap0
	;;-----
	;; Return domain list
	;; 		positions are given in kpc unit
	;;-----

	settings	= self->getheader()
	info 	= self->g_info(snap0)

	n_gal 	= N_ELEMENTS(xc2)
	n_mpi 	= info.ncpu

	xc 	= DOUBLE(xc2)*3.086e21 / info.unit_l
	yc 	= DOUBLE(yc2)*3.086e21 / info.unit_l
	zc 	= DOUBLE(zc2)*3.086e21 / info.unit_l
	rr 	= DOUBLE(rr2)*3.086e21 / info.unit_l

	dom_list 	= LONARR(n_gal, n_mpi)

	ftr_name 	= settings.dir_lib + 'fortran/find_domain.so'
		larr = LONARR(20) & darr = DBLARR(20)
		larr(0) 	= n_gal
		larr(1) 	= n_mpi
		larr(2)		= 1L

		darr(0) 	= 1.d

	void 	= CALL_EXTERNAL(ftr_name, 'find_domain', $
		xc, yc, zc, rr, info.hindex, info.levmax, dom_list, larr, darr)

	cut 	= WHERE(dom_list GE 0L, ncut)

	IF ncut GE 1L THEN BEGIN
		RETURN, dom_list
	ENDIF ELSE BEGIN
		self->errorout, '		g_domain: NO DOMAIN LEFT'
		RETURN, 1L
	ENDELSE
END

FUNCTION veluga::g_part, xc2, yc2, zc2, rr2, snap0, dom_list=dom_list, simout=simout

	;;-----
	;; Read Particle within a sphere
	;; IF dom_list is set, read all ptcls in the argued domain list
	;;	/simout 	- output as the raw simulation unit
	;;-----

	settings	= self->getheader()
	num_thread = self.num_thread
	IF ~KEYWORD_SET(dom_list) THEN BEGIN
		dom_list 	= self->g_domain(xc2, yc2, zc2, rr2, snap0)
		dom_list 	= WHERE(dom_list GE 1L) + 1L
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

	ftr_name 	= settings.dir_lib + 'fortran/jsrd_part_totnum.so'
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

	ftr_name 	= settings.dir_lib + 'fortran/jsrd_part.so'
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

	IF ~KEYWORD_SET(simout) THEN BEGIN
		part.xx 	*= (info.unit_l/3.086e21)	;; [kpc]
		part.yy 	*= (info.unit_l/3.086e21)
		part.zz 	*= (info.unit_l/3.086e21)

		part.vx 	*= info.kms
		part.vy 	*= info.kms
		part.vz 	*= info.kms

		part.mp 	*= (info.unit_m / 1.98892d33)

	ENDIF

	TOC, elapsed_time=elt2

    ;PRINT, elt1, ' - totnum ', elt2, ' - read all'
	RETURN, part
END

FUNCTION veluga::g_cell, xc2, yc2, zc2, rr2, snap0, dom_list=dom_list, simout=simout

	;;-----
	;; Read AMR cells within a sphere
	;; xc2, yc2, zc2, rr2 in kpc unit
	;; IF dom_list is set, read all cells in the argued domain list. If not, read all cells inside R + dX
	;;	/simout 	- output as the raw simulation unit
	;;-----

	settings	= self->getheader()
	num_thread = self.num_thread
	IF ~KEYWORD_SET(dom_list) THEN BEGIN
		dom_list 	= self->g_domain(xc2, yc2, zc2, rr2, snap0)
		dom_list 	= WHERE(dom_list GE 1L) + 1L
	ENDIF

	info 	= self->g_info(snap0)
	ncpu 	= N_ELEMENTS(dom_list)
	dir 	= settings.dir_raw + 'output_' + STRING(snap0, format='(I5.5)')

	xc	= xc2 / info.unit_l * 3.086d21
	yc	= yc2 / info.unit_l * 3.086d21
	zc	= zc2 / info.unit_l * 3.086d21
	rr	= rr2 / info.unit_l * 3.086d21

	IF ~KEYWORD_SET(xr) OR ~KEYWORD_SET(yr) OR ~KEYWORD_SET(zr) THEN BEGIN
		xr	= [-1d,1d] * rr + xc
		yr	= [-1d,1d] * rr + yc
		zr	= [-1d,1d] * rr + zc
	ENDIF

	file_a	= dir + '/amr_' + STRING(snap0,'(I5.5)') + '.out'
	file_h	= dir + '/hydro_' + STRING(snap0,'(I5.5)') + '.out'
	file_i	= dir + '/info_' + STRING(snap0,'(I5.5)') + '.txt'
	
	;;-----
	;; MEMORY ALLOCATE
	;;-----
	ftr_name 	= settings.dir_lib + 'fortran/jsamr2cell_totnum.so'
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

	mesh_xg	= DBLARR(ntot,info.ndim)
	mesh_vx	= DBLARR(ntot,info.ndim)
	mesh_dx	= DBLARR(ntot)
	mesh_hd	= DBLARR(ntot,nvarh)
	mesh_lv	= LONARR(ntot)*0L - 10L

	;;-----
	;; READ CELL
	;;-----
	ftr_name	= settings.dir_lib + 'fortran/jsamr2cell.so'
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

		
		void	= CALL_EXTERNAL(ftr_name, 'jsamr2cell', $
			larr, darr, file_a, file_h, file_i, $
			mg_ind, mesh_xg, mesh_dx, mesh_hd, mesh_lv, dom_list)

	;;-----
	;; POST PROCESSING
	;;-----
	cut	= WHERE(mesh_lv GE 0L,ncell)
	cell 	= self->allocate(ncell, type='cell')

	mesh_xg 	= mesh_xg(cut,*)
	mesh_hd 	= mesh_hd(cut,*)
	mesh_lv 	= mesh_lv(cut)
	mesh_dx 	= mesh_dx(cut)

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

	cell.UE 	= mesh_hd(*,4)/mesh_hd(*,0)/(5.d/3.-1.d) * info.unit_T2 / (1.66d-24) * 1.38049d-23 * 1e-3 ;; [km/s]^2


	IF ~KEYWORD_SET(simout) THEN BEGIN

		cell.xx 	*= (info.unit_l/3.086e21)	;; [kpc]
		cell.yy 	*= (info.unit_l/3.086e21)
		cell.zz 	*= (info.unit_l/3.086e21)

		cell.vx 	*= info.kms					;; [km/s]
		cell.vy 	*= info.kms
		cell.vz 	*= info.kms

		cell.dx 	*= (info.unit_l/3.086e21)	;; [kpc]
		cell.mp 	= mesh_hd(*,0) * info.unit_d * (mesh_dx*info.unit_l)^3.d / 1.98892d33 	;; [Msun]

		cell.den 	*= (info.unit_nH)			;; [/cc]
		cell.temp 	*= (1.d/mesh_hd(*,0) * info.unit_T2)	;; [K/mu]

		;cell.p_thermal 	= mesh_hd(*,4) * info.unit_m / info.unit_l / info.unit_t^2 / 1.3806200e-16
	ENDIF

	RETURN, cell

END


FUNCTION veluga::g_cfrac, xc, yc, zc, rr, snap0, aperture=aperture
	;;-----
	;; Compute contamination fractions
	;;		positions are given in kpc unit
	;;		aperture in R_eff unit
	;;-----

	settings	= self->getheader()
	num_thread 	= self.num_thread
	
	IF ~KEYWORD_SET(aperture) THEN aperture = [1.d]

	n_gal 	= N_ELEMENTS(xc)
	n_aper	= N_ELEMENTS(aperture)

	;;------
	;; Get Domain
	;;------
	dom_list 	= self->g_domain(xc, yc, zc, rr*MAX(aperture), snap0)

	cut 	= WHERE(dom_list GE 1L)
	cut2 	= (ARRAY_INDICES(dom_list, cut))(1,*)
	dom_all	= self->g_unique(cut2) + 1L

	;;-----
	;; Read all ptcls
	;;-----
	part 	= self->g_part(0.d, 0.d, 0.d, 0.d, snap0, dom_list=dom_all)

	dm_ind 	= WHERE(part.family EQ 1L, nn_dm)
	dm_xp 	= part.xp(dm_ind,*)
	dm_mm 	= part.mp(dm_ind)

	;;-----
	;; CFrac computation
	;;-----
	info 	= self->g_info(snap0)

	conf_n 	= DBLARR(n_gal, n_aper)
	conf_m 	= DBLARR(n_gal, n_aper)

	xc0 	= xc * 3.086d21/info.unit_l
	yc0 	= yc * 3.086d21/info.unit_l
	zc0 	= zc * 3.086d21/info.unit_l
	rr0 	= rr * 3.086d21/info.unit_l

	dmp_mass 	= 1./(settings.neff*1.d)^3 * (info.omega_M - info.omega_b) / info.omega_m

	ftr_name 	= settings.dir_lib + 'fortran/get_contam.so'
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
		larr, darr, settings.dir_raw, xc0, yc0, zc0, rr0, $
		dm_xp, dm_mm, conf_n, conf_m, DOUBLE(aperture))
	
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

	ftr_name	= settings.dir_lib + 'fortran/prop_time.so'
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

FUNCTION veluga::g_luminosity, mp, ap, zp, band

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

		FOR j=1L, N_ELEMENTS(ref_met)-1L DO BEGIN
			cut 	= WHERE(met_ind EQ -1L AND met_arr - ref_met(j) LT 0., ncut)
			IF ncut GE 1L THEN met_ind(cut)	= j-1
		ENDFOR

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
		ENDCASE
	ENDFOR

	RETURN, flux
END

FUNCTION veluga::g_mag, lum, band

	
	FOR i=0L, N_ELEMENTS(band)-1L DO BEGIN
		CASE band(i) OF
			'u'		: mag 	= 6.55 - 2.5 * ALOG10(TOTAL(lum))
			'g'		: mag 	= 5.12 - 2.5 * ALOG10(TOTAL(lum))
			'r'		: mag 	= 4.68 - 2.5 * ALOG10(TOTAL(lum))
			'i'		: mag 	= 4.57 - 2.5 * ALOG10(TOTAL(lum))
			'z'		: mag 	= 4.54 - 2.5 * ALOG10(TOTAL(lum))
			'nuv'	: mag 	= 10.18 - 2.5 * ALOG10(TOTAL(lum))
		ENDCASE
	ENDFOR

	RETURN, mag
END

PRO veluga::g_potential, cell, part, $
	p_type=p_type, e_type=e_type, bsize=bsize

	;;-----
	;; Compute potential using all mass components
	;;
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

	npart 	= N_ELEMENTS(part)
	ncell 	= N_ELEMENTS(cell)

	pos 	= DBLARR(npart+ncell,3)
	mass 	= DBLARR(npart+ncell)

	pos(0L:npart-1L,0) 	= part.xx
	pos(0L:npart-1L,1) 	= part.yy
	pos(0L:npart-1L,2) 	= part.zz
	mass(0L:npart-1L) 	= part.mp

	cut_part	= WHERE(part.family NE 1L AND part.family NE 2L, ntracer)
	IF ntracer GE 1L THEN mass(cut_part) = 0.d

	pos(npart:npart+ncell-1L,0) 	= cell.xx
	pos(npart:npart+ncell-1L,1) 	= cell.yy
	pos(npart:npart+ncell-1L,2) 	= cell.zz
	mass(npart:npart+ncell-1L)	 	= cell.mp

	pot 	= DBLARR(npart+ncell)
	force 	= DBLARR(npart+ncell)
	IF N_ELEMENTS(mass) LE bsize THEN bsize = 4L

	STOP
	;;-----
	;; Compute potential
	;;-----
	Gconst 		= 6.67408d-11 		;; m^3 kg^-1 s^-2
	mtokpc 		= (1./3.086d19)
	kgtoMsun 	= (1./1.98892d30)
	Gconst 		*= (mtoKpc / kgtoMsun * 1d-6)	;; (km/s)^2 Kpc/Msun


	ftr_name 	= settings.dir_lib + 'fortran/js_getpt_ft.so'
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

	part.PE 	= pot(0L:npart-1L) 		;; [km/s]^2
	cell.PE 	= pot(npart:npart+ncell-1L)

	;part(cut_part).KE 	= self->g_d3d( part(cut_part).vx, part(cut_part).vy, part(cut_part).vz, []
	RETURN
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
	tbl.metal 	= (10.^metal)*0.02
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
	tbl.metal 	= (10.^metal)*0.02
	tbl.nuv 	= nuv2

	age_arr		= js_unique(tbl.age)
	met_arr		= js_unique(tbl.metal)

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
