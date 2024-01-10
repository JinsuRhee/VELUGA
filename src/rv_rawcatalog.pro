FUNCTION rv_RawCatalog, settings, veluga, runstat, run=run

dir_data	= runstat.dir

;;-----
;; Check procedure set
;;-----
IF run EQ 0L THEN RETURN, PTR_NEW({rv_rawcatalog:-1},/no_copy)
IF run EQ 1L THEN BEGIN
	IF STRLEN(FILE_SEARCH(dir_data + '/rv_io.sav')) GE 5L THEN BEGIN
		RESTORE, dir_data + '/rv_io.sav'
		RETURN, PTR_NEW(output,/no_copy)
	ENDIF ELSE BEGIN
		run	= 2L
	ENDELSE
ENDIF
IF run EQ 2L THEN BEGIN
	;PRINT, '        %%%%% (No previous works are found)'

	;;-----
	;; Find # of columns first, tags, and types
	;;-----
	str	= ' '
	IF settings.horg EQ 'h' THEN OPENR, 10, dir_data + '/halo.dat.properties.0'
	IF settings.horg EQ 'g' THEN OPENR, 10, dir_data + '/galaxy.dat.properties.0'

	READF, 10, str & READF, 10, str & READF, 10, str
	dtag    = STRSPLIT(str, ' ', /extract)
	n_column= N_ELEMENTS(dtag)
	dtype   = REPLICATE('D', n_column)
 	for i=0L, n_column-1L do begin
 	        dum     = STRPOS(dtag(i),'(')
 	        dtag(i) = STRMID(dtag(i),0, dum)
 	        IF dtag(i) EQ 'ID' 		THEN dtype(i) = 'LL'
 	        IF dtag(i) EQ 'ID_mbp'		THEN dtype(i) = 'LL'
 	        IF dtag(i) EQ 'ID_minpot'	THEN dtype(i) = 'LL'
 	        IF dtag(i) EQ 'hostHaloID'	THEN dtype(i) = 'LL'
 	        IF dtag(i) EQ 'numSubStruct'	THEN dtype(i) = 'L'
 	        IF dtag(i) EQ 'npart'		THEN dtype(i) = 'L'
 	        IF dtag(i) EQ 'Structuretype'	THEN dtype(i) = 'L'
 	        IF dtag(i) EQ 'n_gas'		THEN dtype(i) = 'L'
 	        IF dtag(i) EQ 'n_star'		THEN dtype(i) = 'L'
 	        IF dtag(i) EQ 'n_bh'		THEN dtype(i) = 'L'
 	ENDFOR
 	CLOSE, 10

	;;-----
	;; Read Data w/ calling readcol
	;;-----

	IF settings.horg EQ 'h' THEN dum_fname = file_search(dir_data + '/halo.dat.properties.*')
	IF settings.horg EQ 'g' THEN dum_fname = file_search(dir_data + '/galaxy.dat.properties.*')

	dum_nn	= 0L
	FOR i=0L, N_ELEMENTS(dum_fname)-1L DO dum_nn = $
		dum_nn + FILE_LINES(dum_fname(i)) - 3L

	;;;;----- Allocate Memory
	FOR i=0L, N_ELEMENTS(settings.column_list)-1L DO BEGIN
		tmp     = 'arr' + STRTRIM(i,2) + '='
		cut     = WHERE(dtag eq settings.column_list(i), ncut)

		IF ncut EQ 0L THEN PRINT, '        ***** column_list has a wrong argument'
		IF ncut EQ 0L THEN PRINT, '        ***** 	STOP AT rv_RawCatalog.pro'
		IF ncut EQ 0L THEN STOP

		IF dtype(cut) EQ 'LL'THEN tmp = tmp + 'LON64ARR('
		IF dtype(cut) EQ 'L' THEN tmp = tmp + 'LONARR('
		IF dtype(cut) EQ 'D' THEN tmp = tmp + 'DBLARR('
		IF dtype(cut) EQ 'F' THEN tmp = tmp + 'FLTARR('
		tmp     = tmp + 'dum_nn)'
		void    = EXECUTE(tmp)

		IF i EQ 0L THEN tmp2 = 'output = CREATE_STRUCT(dtag(' + $
			STRTRIM(cut(0),2) + '), arr' + STRTRIM(i,2)
		IF i GE 1L THEN tmp2 = tmp2 + ',dtag(' + STRTRIM(Cut(0),2) + $
			'), arr' + STRTRIM(i,2)
	endfor
	tmp2    = tmp2 + ')'
	void    = EXECUTE(tmp2)

	;;;;----- Read
	n1 = 0L &  n2 = -1L
	FOR fi=0L, N_ELEMENTS(dum_fname)-1L DO BEGIN
		fname2  = dum_fname(fi)

		IF FILE_LINES(fname2) EQ 3L THEN CONTINUE
		dum_str = 'READCOL, fname2, '
		FOR i=1L, n_column DO dum_str = dum_str + 'v' + STRTRIM(i,2) + ', '
		dum_str = dum_str + 'format="'
		FOR i=0L, n_column-2L DO dum_str = dum_str + STRTRIM(dtype(i),2) + ', '
		dum_str = dum_str + STRTRIM(dtype(n_column-1L),2) + '", '
		dum_str = dum_str + 'numline=file_lines(fname2), skipline=3, /silent'
		void    = EXECUTE(dum_str)

		;;;;;;----- Adjust the line numbers
		n1	= n2 + 1L
		n2      = n2 + N_ELEMENTS(v1)

		;;;;;;----- Extract the requested columns
		n_match	= 0L
		n_nan	= 0L
		FOR i=0L, n_column-1L DO BEGIN
			cut     = WHERE(settings.column_list EQ dtag(i), ncut)
			IF ncut EQ 0L THEN CONTINUE
			n_match ++
			dum_str = 'output.' + STRTRIM(dtag(i),2) + $
				'(' + STRTRIM(n1,2) + ':' + STRTRIM(n2,2) + ') = ' + $
				'v' + STRTRIM(i+1,2)
			void    = EXECUTE(dum_str)
		ENDFOR

		IF N_ELEMENTS(settings.column_list) NE n_match THEN $
			PRINT, '        ***** There is a wrong typed one in column_list'
	ENDFOR

	IF settings.pp_saveprocess EQ 1L THEN SAVE, filename=dir_data + '/rv_io.sav', output
	RETURN, PTR_NEW(output,/no_copy)
ENDIF
END
