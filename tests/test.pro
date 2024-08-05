PRO test

	v 	= obj_new('veluga', 'vrheader_ref.txt', num_thread=10L)
	c 	= v->g_cell(811L, 0.d, 0.d, 0.d, 0.d, dom_list=LINDGEN(100)+1L)

	
	STOP
END