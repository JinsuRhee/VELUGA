PRO halo_visualcheck, v, horg=horg, snap=snap, outdir=outdir, mrange=mrange, rfact=rfact

	IF ~KEYWORD_SET(horg) THEN BEGIN
		horg = ''
		READ, horg, prompt='input halo/galaxy identifier) halo(h) / galaxy(g): '
	ENDIF
	IF ~KEYWORD_SET(snap) THEN BEGIN
		snap = 0L
		READ, snap, prompt='input snapshot number: '
	ENDIF
	IF ~KEYWORD_SET(outdir) THEN BEGIN
		outdir = ''
		READ, outdir, prompt='path for saving outputs: '
	ENDIF
	IF ~KEYWORD_SET(mrange) THEN BEGIN
		mlow	= 0.d
		mhigh	= 0.d
		READ, mlow, prompt='lower mass cut: '
		READ, mhigh, prompt='upper mass cut: '
		mrange	= [mlow, mhigh]
	ENDIF
	IF ~KEYWORD_SET(rfact) THEN BEGIN
		rfact	= 0.d
		READ, rfact, prompt='size of the box in radius unit (if a negative value given, the fixed size is used: '
	ENDIF

	;;-----
	;; LOAD GAL
	;;-----
	halo	= v->r_gal(snap, -1L, horg=horg)
	IF horg EQ 'h' THEN mass = halo.mass_200crit ELSE mass = halo.mass_tot
	cut	= WHERE(mass GE mrange(0) AND mass LT mrange(1), ncut)
	IF ncut EQ 0L THEN BEGIN
		PRINT, 'no sample left with this mass range : ', mrange
		STOP
	ENDIF
	halo0	= halo
	halo	= halo(cut)

	;;-----
	;; DRAW	
	;;-----
	color	= ['red', 'light coral', 'orange', 'yellow', 'green', 'cyan', $
		'sea green', 'dodger blue', 'blue', 'navy', 'purple', 'brown']
	FOR i=0L, N_ELEMENTS(halo)-1L DO BEGIN
;IF i NE N_ELEMENTS(halo)-1L THEN CONTINUE
		tmp	= halo(i)
		IF rfact GT 0. THEN BEGIN
			IF horg EQ 'h' THEN rr = tmp.r_200crit*rfact ELSE rr = tmp.r_halfmass*rfact
		ENDIF ELSE BEGIN
			rr = ABS(rfact)
		ENDELSE
		;; set box
		xr	= [-1.d,1.d]*rr + tmp.xc
		yr	= [-1.d,1.d]*rr + tmp.yc
		zr	= [-1.d,1.d]*rr + tmp.zc

		;; find subhalos
		inbox	= v->g_boundind(xx=halo0.xc, yy=halo0.yc, zz=halo0.zc, xr=xr, yr=yr, zr=zr)
		If inbox.n EQ 0L THEN CONTINUE
		sat	= halo0(inbox.ind)
		remhost	= WHERE(sat.ID NE tmp.id)
		sat	= sat(remhost)
		sat	= sat( REVERSE(SORT(sat.mass_tot)) )

		;; gather all particles
		p_all	= v->g_part(snap, tmp.xc, tmp.yc, tmp.zc, rr*2.d, /memeff)
		IF horg EQ 'h' THEN famtype = 1L ELSE famtype = 2L
		cut	= WHERE(*p_all.family EQ famtype)
		p_all	= v->g_extract(p_all, cut)


		IF horg EQ 'h' THEN BEGIN ;; remove poor resolution DM
			info 	= v->g_info(snap)
			settings	= v->getheader()
			dmp_mass        = 1.d/(settings.neff*1.d)^3 * (info.omega_M - info.omega_b)/info.omega_M
			dmp_mass	*= (info.unit_m / info.cgs.m_sun)
			cut	= WHERE(*p_all.mp LT 1.1*dmp_mass)
			p_all	= v->g_extract(p_all, cut)
		ENDIF

		inbox	= v->g_boundind(xx=*p_all.xx, yy=*p_all.yy, zz=*p_all.zz, xr=xr, yr=yr, zr=zr)
		p_all0	= p_all
		p_all	= v->g_extract(p_all, inbox.ind)

		;; load host particles
		pid	= v->r_pid(snap, tmp.id, horg=horg)
		idmat	= v->g_indmatch(*p_all.id, pid)
		cut	= WHERE(idmat.x GE 0L)
		p_host	= v->g_extract(p_all, cut)

		;; load sats
		pid	= []
		refsat	= PTRARR(12)
		nref	= 0L
		FOR j=0L, N_ELEMENTS(sat)-1L DO BEGIN
			pid0	= v->r_pid(snap, sat(j).id, horg=horg)
			pid	= [pid, pid0]

			IF j LE N_ELEMENTS(refsat)-1L THEN BEGIN
				idmat	= v->g_indmatch(*p_all0.id, pid0)
				cut	= WHERE(idmat.x GE 0L)
				refsat(j)	= PTR_NEW(v->g_extract(p_all0, cut))
				nref ++
			ENDIF
		ENDFOR
		idmat	= v->g_indmatch(*p_all.id, pid)
		cut	= WHERE(idmat.x GE 0L)
		p_satall= v->g_extract(p_all, cut)

		;; Get Denstiy
		d_all	= (v->d_2dmap(*p_all.xx, *p_all.yy, xr=xr, yr=yr, n_pix=1000L, mode=-1L, kernel=1L, zz=*p_all.mp)).z
		d_host	= (v->d_2dmap(*p_host.xx, *p_host.yy, xr=xr, yr=yr, n_pix=1000L, mode=-1L, kernel=1L, zz=*p_host.mp)).z
		d_satall= (v->d_2dmap(*p_satall.xx, *p_satall.yy, xr=xr, yr=yr, n_pix=1000L, mode=-1L, kernel=1L, zz=*p_satall.mp)).z
		d_rem	= (v->d_2dmap($
			[*p_all.xx, *p_host.xx, *p_satall.xx], $
			[*p_all.yy, *p_host.yy, *p_satall.yy], $
			zz=[*p_all.mp, (*p_host.mp)*(-1.d), (*p_satall.mp)*(-1.d)], $
			xr=xr, yr=yr, n_pix=1000L, mode=-1L, kernel=1L)).z

		d_sat	= DBLARR(1000L, 1000L, N_ELEMENTS(refsat))
		FOR j=0L, N_ELEMENTS(refsat)-1L DO BEGIN
			IF j EQ nref THEN BREAK

			IF horg EQ 'h' THEN rr2 = sat(j).r_200crit*1.5 ELSE rr2 = sat(j).r_halfmass*1.5

			xr2	= [-1.d, 1.d]*rr2 + sat(j).xc
			yr2	= [-1.d, 1.d]*rr2 + sat(j).yc
			den0	= v->d_2dmap( *(*refsat(j)).xx, *(*refsat(j)).yy, zz=*(*refsat(j)).mp, xr=xr2, yr=yr2, n_pix=1000L, mode=-1L, kernel=1L)
			d_sat(*,*,j)	= den0.z
		ENDFOR

		;; Get Image
		i_all	= BYTE(255.d * v->d_minmax(d_all, 1e-5*MAX(d_all), MAX(d_all), stype='log2') )
		i_host	= BYTE(255.d * v->d_minmax(d_host, 1e-5*MAX(d_host), MAX(d_host), stype='log2') )
		i_satall	= BYTE(255.d * v->d_minmax(d_satall, 1e-5*MAX(d_satall), MAX(d_satall), stype='log2') )
		i_rem	= BYTE(255.d * v->d_minmax(d_rem, 1e-5*MAX(d_rem), MAX(d_rem), stype='log2') )

		i_sat	= BYTARR(1000L, 1000L, N_ELEMENTS(refsat))
		FOR j=0L, N_ELEMENTS(refsat)-1L DO BEGIN
			i_sat(*,*,j) = $
				BYTE(255.d * v->d_minmax(d_sat(*,*,j), 1e-5 * MAX(d_sat(*,*,j)), MAX(d_sat(*,*,j)), stype='log2') )
		ENDFOR

		;; Draw
		iname	= outdir + '/I_' + STRING(tmp.ID, format='(I6.6)') + '.eps'
		cgPS_open, iname, /encapsulated
		
		cgDisplay, 1200., 700.

		cgImage, i_all  , position=[0.,   2./3.5, 0.25, 1.], /noerase
		cgImage, i_host , position=[0.25, 2./3.5, 0.5, 1.], /noerase
		cgImage, i_satall,position=[0.5,  2./3.5, 0.75, 1.], /noerase
		cgImage, i_rem,   position=[0.75, 2./3.5, 1., 1.], /noerase

		pos0	= [0., 1./3.5, 1./6, 2./3.5]
		FOR j=0L, N_ELEMENTS(refsat)-1L DO BEGIN
			cgImage, i_sat(*,*,j), position=pos0, /noerase
			cgText, pos0(0)+0.01, pos0(3)-0.03, 'Sat ' + STRING(j,format='(I2)'), charsize=0.8, charthick=3.0, color=color(j), /normal
			pos0	+= [1./6, 0.d, 1./6, 0.d]
			IF j EQ 5L THEN BEGIN
				pos0	= [0., 0., 1./6, 1./3.5]
			ENDIF
		ENDFOR

		cgPlot, 0, 0, /nodata, /noerase, position=[0.5, 2./3.5, 0.75, 1.], xr=xr, yr=yr, xstyle=4, ystyle=4
		FOR j=0L, N_ELEMENTS(refsat)-1L DO BEGIN
			cgOplot, sat(j).xc, sat(j).yc, psym=16, symsize=0.25, color=color(j)
		ENDFOR

		!p.font = -1
		dtx	= 0.02
		dty	= -0.04
		cgText, 0.00 + dtx, 1. + dty, 'All', charsize=1.2,       color='red', charthick=3.5, /normal
		cgText, 0.25 + dtx, 1. + dty, 'Host only', charsize=1.2, color='red', charthick=3.5, /normal
		cgText, 0.50 + dtx, 1. + dty, 'Sats only', charsize=1.2, color='red', charthick=3.5, /normal
		cgText, 0.75 + dtx, 1. + dty, 'Residual', charsize=1.2,  color='red', charthick=3.5, /normal

		;cgText, 0.00 + dtx, 0.0 + dty, 'Sat 1', charsize=1.2, color=color(0), charthick=3.5, /normal
		;cgText, 0.25 + dtx, 0.0 + dty, 'Sat 2', charsize=1.2, color=color(1), charthick=3.5, /normal
		;cgText, 0.50 + dtx, 0.0 + dty, 'Sat 3', charsize=1.2, color=color(2), charthick=3.5, /normal
		;cgText, 0.75 + dtx, 0.0 + dty, 'Sat 4', charsize=1.2, color=color(3), charthick=3.5, /normal
		cgPS_close


		v->free, p_all
		v->free, p_all0
		v->free, p_host
		v->free, p_satall
		PTR_FREE, refsat

;		STOP
	ENDFOR
	STOP
END
