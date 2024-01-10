Pro simple_write_hdf5, data, name, fid, skip=skip

	IF KEYWORD_SET(skip) THEN RETURN

	if max(size(data, /dimensions)) eq 0L then data = [data]

	datatype_id     = h5t_idl_create(data)
	dataspace_id    = h5s_create_simple(size(data, /dimensions))
	dataset_id      = h5d_create(fid, $
		name, datatype_id, dataspace_id)

	h5d_write, dataset_id, data
	h5d_close, dataset_id
	h5s_close, dataspace_id
	h5t_close, datatype_id

End
