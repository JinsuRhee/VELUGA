#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include "hdf5.h"

typedef struct {
   unsigned short slen;         /* length of the string         */
   short stype;                 /* Type of string               */
   char *s;                     /* Pointer to chararcter array  */
} STRING;

#define STR_LEN(__str)    ((long)(__str)->slen)

void get_data(hid_t gid, double* d_array, int* l_array, char* base, int* header, int n_gal, int d_nn, bool isdouble, bool isgprop){


	char name[256];
	if(isgprop) snprintf(name, sizeof(name), "G_%s", base);
	else snprintf(name, sizeof(name), "%s", base);
	

	hid_t dset_id, space_id;
	herr_t status;
	dset_id	= H5Dopen(gid, name, H5P_DEFAULT);
	space_id = H5Dget_space(dset_id);

	hsize_t	dims[3] = {0};
	H5Sget_simple_extent_dims(space_id, dims, NULL);

	if(d_nn != dims[0]){
		printf("Wrong Mapping: %s // d_nn = %d and dims = %d \n ", base, d_nn, dims[0]);
	}

	if(isdouble){
		double *buf = (double *)malloc(sizeof(double) * dims[0]);
		status	= H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

		for(int i=0; i<dims[0]; i++){
			//d_array[(*header) + i*n_gal] = buf[i];
			d_array[(*header)] = buf[i];
			//printf("INDEX = %d \n", (*header) + i*n_gal);
			//d_array[*header] = buf[i];
			
			(*header) ++;
		};

		
		

		free(buf);
	}else{
		int *buf = (int *)malloc(sizeof(int) * dims[0]);
		status	= H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

		for(int i=0; i<dims[0]; i++){
			//l_array[(*header) + i*n_gal] = buf[i];
			l_array[(*header)] = buf[i];
			//printf("INDEX = %d \n", (*header) + i*n_gal);
			//l_array[*header] = buf[i];
			
			(*header) ++;
		};

		
		
		free(buf);
	}

	
	H5Sclose(space_id);
	H5Dclose(dset_id);

}


void read_cat(int argc, void *argv[])
{

  STRING *fname, *flux_list, *gprop_map;
  int *larr;
  double *darr;

  int *snap0, *id0, *id, *l_array, *gprop_type, *gprop_nn, *gprop_tag;
  double *d_array;


  //-----
  // Receive data
  //-----

  //----- INFO
  larr		= (int *) argv[0];
  darr		= (double *) argv[1];
  fname		= (STRING *) argv[2];

  snap0		= (int *) argv[3];
  id0		= (int *) argv[4];
  id		= (int *) argv[5];
  d_array	= (double *) argv[6];
  l_array	= (int *) argv[7];
  gprop_tag	= (int *) argv[8];
  gprop_type	= (int *) argv[9];
  gprop_nn	= (int *) argv[10];
  gprop_map	= (STRING *) argv[11];
  flux_list		= (STRING *) argv[12];

  // elements number
  int n_gal		= larr[0];
  int nd_double		= larr[1];
  int nd_long	 	= larr[2];
  int n_prop		= larr[3];
  int n_flux		= larr[4];
  int n_all 		= larr[5];
  int n_propall		= larr[6];
  int n_aper		= larr[7];

  //-----
  // Set pointer
  //-----
  int wpoint[n_prop];
  int wpoint2, wpoint3;
  for(int i=0; i<n_prop; i++){
	  wpoint[i] = 0;
  }

  int int_p=0;
  int dbl_p=0;
  for(int i=0; i<n_prop; i++){

	  if(gprop_type[i] == 1){
		  wpoint[i]	= int_p;
		  int_p += gprop_nn[i]*n_gal;
	  }
	  if(gprop_type[i] == 2){
		  wpoint[i]	= dbl_p;
		  dbl_p += gprop_nn[i]*n_gal;
  	}

  	//printf(" ind : %s %d %d %d \n", gprop_map[gprop_tag[i]].s, wpoint[i], int_p, dbl_p);
  }

  //-----
  // Mapping
  //-----

  //-----
  // Open HDF5
  //-----
  hid_t	file_id, dataset_id, dataspace_id; // hdf5 identifier

  file_id = H5Fopen(fname->s, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0) return;

  //-----
  // READ
  //-----
  char idbase[64];
  char idbase2[64];
  char formatted[64];
  char propbase[64];
  //char propbase2[64];


  for(int i=0; i<n_all; i++){
	  if(*id0>0 && *id0 != id[i]) continue;

	  // Open Group
	  snprintf(idbase, sizeof(idbase), "%s_%06d", "/ID", id[i]);
	  //sprintf(idbase, "/ID_")
	  //sprintf(formatted, "%06d", id[i]);
	  //strcat(idbase, formatted);
	  hid_t group_id = H5Gopen(file_id, idbase, H5P_DEFAULT);

	  snprintf(idbase2, sizeof(idbase), "%s_%06d%s", "/ID", id[i], "/G_Prop");
	  //strcat(idbase2, formatted);
	  //strcat(idbase2, "/G_Prop");
	  hid_t group_id2 = H5Gopen(file_id, idbase2, H5P_DEFAULT);


	  for(int j=0; j<n_prop; j++){

		  //propbase	= gprop_map[gprop_tag[j]].s;
	  	  //printf(" !!!  %d / %d \n", j, wpoint[i]);

		  if(gprop_tag[j] == 32 || gprop_tag[j] == 33){
			  //propbase	= 

			  
			  //printf(" ?? %d : %d \n", j, wpoint[j]);
			  for(int k=0; k<n_flux; k++){
			  	  wpoint2	= wpoint[j] + n_gal*n_aper*k;


				  //char propname[64];
				  snprintf(propbase, sizeof(propbase), "%s_%s", gprop_map[gprop_tag[j]].s, flux_list[k].s);
				  get_data(group_id2, d_array, l_array, propbase, &wpoint2, n_gal, gprop_nn[j]/n_flux, true, true);


				  //printf(" @@@@@ %d / %d \n", k, wpoint2);

				  if(k==0) wpoint3	= wpoint2;
			  }

			  wpoint[j]	= wpoint3;
			  
		  }else if(gprop_tag[j] == 36 || gprop_tag[j] == 37){
			  //propbase	= gprop_map[gprop_tag[j]].s;
			  snprintf(propbase, sizeof(propbase), "%s", gprop_map[gprop_tag[j]].s);
			  get_data(group_id, d_array, l_array, propbase, &wpoint[j], n_gal, gprop_nn[j], false, false);
		  }else{
			  //propbase	= gprop_map[gprop_tag[j]].s;
			  snprintf(propbase, sizeof(propbase), "%s", gprop_map[gprop_tag[j]].s);
			  if(gprop_type[j] == 1){
				  get_data(group_id2, d_array, l_array, propbase, &wpoint[j], n_gal, gprop_nn[j], false, true);
			  }else if(gprop_type[j] == 2){
				  get_data(group_id2, d_array, l_array, propbase, &wpoint[j], n_gal, gprop_nn[j], true, true);
			  }
		  }

		  //printf(" !!! %d / %s / %d \n", i, gprop_map[gprop_tag[j]].s, wpoint[j]);

	  }
	  H5Gclose(group_id);
	  H5Gclose(group_id2);
  }
  //-----
  // Close
  //-----
  H5Fclose(file_id);

}


