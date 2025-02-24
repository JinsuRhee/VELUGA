PRO run_ctree_integrity, vrheader, s0, s1, horg=horg

	IF ~KEYWORD_SET(horg) THEN BEGIN
		horg = 'g'
	ENDIF

	veluga	= OBJ_NEW('veluga', vrheader, num_thread=1L)

	settings 	= veluga.getheader()

	IF horg EQ 'g' THEN BEGIN
		tname 	= settings.dir_catalog + 'Galaxy/tree/tfout/ctree.sav'
	ENDIF ELSE BEGIN
		tname  	= settings.dir_catalog + 'Halo/tree/tfout/ctree.sav'
	ENDELSE

	RESTORE, tname


	nb	= N_ELEMENTS(complete_tree)

	FOR i=s0, s1 DO BEGIN

		g	= veluga->r_gal(i, -1L, horg=horg, Gprop=['ID'])
		IF TYPENAME(g) EQ 'LONG' THEN CONTINUE

		FOR j=0L, N_ELEMENTS(g)-1L DO BEGIN
			i0 	= g(j).id
			s0	= g(j).snapnum

			tk	= s0+i0*tree_key(0)
			tk	= tree_key(tk)
			IF tk LT 0L THEN CONTINUE
			ct0	= *complete_tree(tk)


			;; Key Test

			klist 	= ct0.snap + ct0.ID*tree_key(0)
			klist	= tree_key(klist)

			void 	= WHERE(klist NE tk, nn)

			IF nn GE 1L THEN STOP

			;; Refer from others?
			void 	= WHERE(tree_key EQ tk, nn2)

			IF nn2 NE N_ELEMENTS(klist) THEN STOP
		ENDFOR
		PRINT, i, ' / ', s1
	ENDFOR
END
