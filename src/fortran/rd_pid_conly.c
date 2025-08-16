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

void save_cat(int argc, void *argv[])
{

  STRING *fname, *flux_list;
  int *larr;
  double *darr;

  double *sfr_r, *sfr_t, *mag_r, *conf_r;

  //-----
  // Receive data
  //-----

  //----- INFO
  larr		= (int *) argv[0];
  darr		= (double *) argv[1];
  fname		= (STRING *) argv[2];

  sfr_r			= (double *) argv[3];
  sfr_t			= (double *) argv[4];
  mag_r			= (double *) argv[5];
  conf_r		= (double *) argv[6];

  flux_list		= (STRING *) argv[7];

  //----- Bulk
  double *b_mass_tot, *b_r_halfmass, *b_mvir, *b_rvir, *b_m200, *b_r200, *b_cm, *b_cn;
  int *b_id;
  b_mass_tot		= (double *) argv[8];
  b_r_halfmass		= (double *) argv[9];
  b_mvir		= (double *) argv[10];
  b_rvir		= (double *) argv[11];
  b_m200		= (double *) argv[12];
  b_r200		= (double *) argv[13];
  b_id			= (int *) argv[14];
  b_cm			= (double *) argv[15];
  b_cn			= (double *) argv[16];

  //----- Gal/Particle ID
  long long *p_bid, *p_uid, *p_id;
  p_bid			= (long long *) argv[17];
  p_uid			= (long long *) argv[18];
  p_id			= (long long *) argv[19];


  //----- Catalog Data
  int *cat_long;
  double *cat_dbl;
  STRING *clist_long, *clist_dbl;

  cat_long		= (int *) argv[20];
  clist_long		= (STRING *) argv[21];
  cat_dbl		= (double *) argv[22];
  clist_dbl		= (STRING *) argv[23];

  //----- Bulk Properties
  double *b_sfr, *b_mag, *b_sb, *b_cfm, *b_cfn;
  int *b_isclump, *b_domlist;
  b_sfr			= (double *) argv[24];
  b_mag			= (double *) argv[25];
  b_sb			= (double *) argv[26];
  b_cfm			= (double *) argv[27];
  b_cfn			= (double *) argv[28];
  b_isclump		= (int *) argv[29];
  b_domlist		= (int *) argv[30];

  int *dum		= (int *) argv[30];
  // elements number
  int n_gal		= larr[0];
  int n_mpi		= larr[1];
  int n_sfr_aper 	= larr[10];
  int n_mag_aper	= larr[11];
  int n_conf_aper	= larr[12];
  int n_flux		= larr[13];
  int n_bind		= larr[14];
  int n_catl		= larr[15];
  int n_catd		= larr[16];
  int n_bsfr		= larr[17];
  int n_bmag		= larr[18];
  int n_bsb		= larr[19];
  int n_isclump		= larr[20];
  int n_ndom		= larr[21];

  double aexp		= darr[0];
  //-----
  // Open HDF5
  //-----
  hid_t	file_id, dataset_id, dataspace_id; // hdf5 identifier

  file_id = H5Fcreate(fname->s, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  //-----
  // SAVE DATA
  //-----
  
  // Some info
  save_cat_write_d(file_id, sfr_r,  n_sfr_aper,  "SFR_R");
  save_cat_write_d(file_id, sfr_t,  n_sfr_aper,  "SFR_T");
  save_cat_write_d(file_id, mag_r,  n_mag_aper,  "MAG_R");
  save_cat_write_d(file_id, conf_r, n_conf_aper, "CONF_R");

  char *flist[n_flux];
  for(int i=0; i<n_flux; i++){flist[i] = flux_list[i].s;}
  save_cat_write_s(file_id, flist, n_flux, "Flux_List");

  //for(int i=0; i<n_flux; i++){printf("%s\n", flist[i]); fflush(stdout);}

  // Bulk Properties
  save_cat_write_d(file_id, b_mass_tot,  n_gal,  "Mass_tot");
  save_cat_write_d(file_id, b_r_halfmass,  n_gal,  "R_HalfMass");
  save_cat_write_d(file_id, b_mvir,  n_gal,  "Mvir");
  save_cat_write_d(file_id, b_rvir,  n_gal,  "Rvir");
  save_cat_write_d(file_id, b_m200,  n_gal,  "Mass_200crit");
  save_cat_write_d(file_id, b_r200,  n_gal,  "R_200crit");
  save_cat_write_i(file_id, b_id,  n_gal,  "ID");
   
  save_cat_write_d2(file_id, b_cm,  n_gal, n_conf_aper, "CONF_M");
  save_cat_write_d2(file_id, b_cn,  n_gal, n_conf_aper, "CONF_N");
 
  // Individual galaxy
  for(int i=0; i<n_gal; i++){
	  //printf("gal = %d / %d\n", i, n_gal-1);
	  //fflush(stdout);

	  // Open Group
	  char idbase[64] = "/ID_";
	  char idbase_p[64] = "/ID_";
	  char idbase_g[64] = "/ID_";
	  char formatted[16];

    	  sprintf(formatted, "%06d", b_id[i]);
    	  strcat(idbase, formatted);

    	  sprintf(formatted, "%06d", b_id[i]);
    	  strcat(idbase_p, formatted);
	  strcat(idbase_p, "/P_Prop");

    	  sprintf(formatted, "%06d", b_id[i]);
    	  strcat(idbase_g, formatted);
	  strcat(idbase_g, "/G_Prop");

	  hid_t group_id = H5Gcreate(file_id, idbase, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  hid_t group_id_p = H5Gcreate(file_id, idbase_p, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  hid_t group_id_g = H5Gcreate(file_id, idbase_g, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  

	  // particle id
	  long long ib[2];
	  long long iu[2]; 
	  long long *ptcl_id;

	  if(n_bind > 1){
	        ib[0] = p_bid[i];
	        ib[1] = p_bid[i+n_gal];
	        iu[0] = p_uid[i];
	        iu[1] = p_uid[i+n_gal];

	        long long numptcl = (iu[1] - iu[0] + 1) + (ib[1] - ib[0] + 1);

	        ptcl_id = (long long *)malloc(sizeof(long long) * numptcl);
	        if (ptcl_id == NULL) {
                	fprintf(stderr, "Memory allocation failed\n");
                	return;
    	        }

	        for(long long j=0; j<ib[1] - ib[0] + 1; j++){
	        	ptcl_id[j] = p_id[ib[0] + j];
	        }

	        for(long long j=0; j<iu[1] - iu[0] + 1; j++){
	        	ptcl_id[j + ib[1] - ib[0]+1] = p_id[iu[0]+j];
	        }

	        long long count = 0;
    	        for (long long j = 0; j < numptcl; j++) {
                	if (ptcl_id[j] > -922337203685477580LL) {
            			count++;
	        	}
	        }

	        //
    	        long long *f_ptcl_id = (long long *)malloc(sizeof(long long) * count);
    	        if (f_ptcl_id == NULL) {
                	fprintf(stderr, "Memory allocation failed for filtered array\n");
                	free(ptcl_id);
                	return;
    	        }

    	        //
    	        long long dumn = 0;
    	        for (long long j = 0; j < numptcl; j++) {
                	if (ptcl_id[j] > -922337203685477580LL) {
            			f_ptcl_id[dumn++] = ptcl_id[j];
                	}
    	        }


	        //write
  		save_cat_write_i64(group_id_p, f_ptcl_id,  count,  "P_ID");

	        free(ptcl_id);
	        free(f_ptcl_id);
	        ptcl_id = NULL;
	        f_ptcl_id = NULL;
	  }

	  //printf("gal pt = %d\n", i);
	  //fflush(stdout);
	  // catalog data (long)
	  int dummy_n = 1;
	  int *dummy_int;
	  for(int j=0; j<n_catl; j++){
	          char catbase[64] = "G_";
	          save_cat_write_i(group_id_g, &cat_long[n_gal*j + i], dummy_n, strcat(catbase,clist_long[j].s));
	  }
	  // catalog data (double)
	  for(int j=0; j<n_catd; j++){
	          char catbase[64] = "G_";
	          save_cat_write_d(group_id_g, &cat_dbl[n_gal*j + i], dummy_n, strcat(catbase, clist_dbl[j].s));
	  }

	  //printf("gal cat = %d\n", i);
	  //fflush(stdout);
	  //// Bulk - SFR
	  double *dummy_sfr;
	  dummy_n = 1;
	  if(n_bsfr >= 2){
	        dummy_sfr = (double *)malloc(sizeof(double) * n_sfr_aper);
	        for(int j=0; j<n_sfr_aper; j++){
	        	dummy_sfr[j] = b_sfr[n_sfr_aper*i+j];
	        }
	        save_cat_write_d(group_id_g, dummy_sfr, n_sfr_aper, "G_SFR");

	        free(dummy_sfr);
	        dummy_sfr = NULL;
	  }else {
	        dummy_sfr = (double *)malloc(sizeof(double) * dummy_n);
	        dummy_sfr[0] = -1.;
	        save_cat_write_d(group_id_g, dummy_sfr, dummy_n, "G_SFR");

	        free(dummy_sfr);
	        dummy_sfr = NULL;
	  }

	  //printf("gal sfr = %d\n", i);
	  //fflush(stdout);
	  //// Bulk - MAG
	  double *dummy_mag;
	  dummy_n = 1;
	  for(int k=0; k<n_flux; k++){
	          char magbase[64] = "G_ABmag_";
	          if(n_bmag >= 2){
	        	dummy_mag = (double *)malloc(sizeof(double) * n_mag_aper);
	        	for(int j=0; j<n_mag_aper; j++){
	        		dummy_mag[j] = b_mag[n_flux*n_mag_aper*i + n_mag_aper*k + j];
	        	}
	        	save_cat_write_d(group_id_g, dummy_mag, n_mag_aper, strcat(magbase,flux_list[k].s));

	        	free(dummy_mag);
	        	dummy_mag = NULL;
	          }else{
	          	dummy_mag = (double *)malloc(sizeof(double) * dummy_n);
	        	dummy_mag[0] = -1.;
	        	save_cat_write_d(group_id_g, dummy_mag, dummy_n, strcat(magbase,flux_list[k].s));

	        	free(dummy_mag);
	        	dummy_mag = NULL;

	          }
	  }

	  //printf("gal mag = %d\n", i);
	  //fflush(stdout);
	  //// Bulk - SB
	  double *dummy_sb;
	  dummy_n = 1;
	  for(int k=0; k<n_flux; k++){
	          char magbase[64] = "G_SB_";
	          if(n_bsb >= 2){
	          	dummy_sb = (double *)malloc(sizeof(double) * n_mag_aper);
	        	for(int j=0; j<n_mag_aper; j++){
	        		dummy_sb[j] = b_sb[n_flux*n_mag_aper*i + n_mag_aper*k + j];
	        	}
	        	save_cat_write_d(group_id_g, dummy_sb, n_mag_aper, strcat(magbase,flux_list[k].s));

	        	free(dummy_sb);
	        	dummy_sb = NULL;
	          } else{
	          	dummy_sb = (double *)malloc(sizeof(double) * dummy_n);
	        	dummy_sb[0] = -1.;
	        	save_cat_write_d(group_id_g, dummy_sb, dummy_n, strcat(magbase,flux_list[k].s));

	        	free(dummy_sb);
	        	dummy_sb = NULL;
	          }
	  }


	  //printf("gal sb = %d\n", i);
	  //fflush(stdout);
	  //// Bulk - Cfrac
	  double *dummy_cfm, *dummy_cfn;
	  dummy_cfm = (double *)malloc(sizeof(double) * n_conf_aper);
	  dummy_cfn = (double *)malloc(sizeof(double) * n_conf_aper);
	  for(int j=0; j<n_conf_aper; j++){
	          dummy_cfm[j] = b_cfm[n_conf_aper*i + j];
	          dummy_cfn[j] = b_cfn[n_conf_aper*i + j];
	  }

	  save_cat_write_d(group_id_g, dummy_cfm, n_conf_aper, "G_ConFrac_M");
	  save_cat_write_d(group_id_g, dummy_cfn, n_conf_aper, "G_ConFrac_N");

	  free(dummy_cfm);
	  free(dummy_cfn);
	  dummy_cfm = NULL;
	  dummy_cfn = NULL;

	  //printf("gal cf = %d\n", i);
	  //fflush(stdout);
	  //// Bulk - isClump
	  dummy_n = 1;
	  if(n_isclump >= 2){
	        save_cat_write_i(group_id, &b_isclump[i], dummy_n, "isclump");
	  }else {
	        save_cat_write_i(group_id, &b_isclump[0], dummy_n, "isclump");

	        free(dummy_mag);
	        dummy_mag = NULL;
	  }

	  //printf("gal ic = %d\n", i);
	  //fflush(stdout);
	  //// Aexp
	  save_cat_write_d(group_id, &aexp, dummy_n, "Aexp");

	  //printf("gal aexp = %d\n", i);
	  //fflush(stdout);
	  //// DOM LIST
	  int *dummy_dom;
	  dummy_n = 1;

	  if(n_ndom >= 2){
	          dummy_dom = (int *)malloc(sizeof(int) * n_mpi);
	          for(int j=0; j<n_mpi; j++){
	        	  dummy_dom[j]	= b_domlist[n_mpi*i + j];
	          }
	          save_cat_write_i(group_id, dummy_dom, n_mpi, "Domain_List");

	          free(dummy_dom);
	          dummy_dom = NULL;
	  }else{
	          save_cat_write_i(group_id, &b_domlist[0], dummy_n, "Domain_List");
	  }
  

	  // confrac value and ptcl id
	  H5Gclose(group_id);
	  H5Gclose(group_id_g);
	  H5Gclose(group_id_p);

  } 

  //-----
  // Close HDF5
  //-----
  H5Fclose(file_id);
}

//for double
void save_cat_write_d(hid_t fid, double *data, int ndata, char *fieldname){

	hid_t dspace_id, dset_id;

	//-----
	// Data SPACE
	//-----
	hsize_t dims[1] = {ndata};
	dspace_id = H5Screate_simple(1, dims, NULL);
	
	//-----
	// Set Space
	//-----
	//printf(fieldname);
	dset_id = H5Dcreate(fid, fieldname, H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//-----
	// SAVE
	//-----
	herr_t status;
	status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	//printf(status);
	//-----
	// Close
	//-----
	H5Dclose(dset_id);
	H5Sclose(dspace_id);
}
//for double 2D
void save_cat_write_d2(hid_t fid, double *data, int ndata, int ndata2, char *fieldname){

	hid_t dspace_id, dset_id;

	//-----
	// Data SPACE
	//-----
	hsize_t dims[2] = {ndata, ndata2};
	dspace_id = H5Screate_simple(2, dims, NULL);
	
	//-----
	// Set Space
	//-----
	//printf(fieldname);
	dset_id = H5Dcreate(fid, fieldname, H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//-----
	// SAVE
	//-----
	herr_t status;
	status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	//printf(status);
	//-----
	// Close
	//-----
	H5Dclose(dset_id);
	H5Sclose(dspace_id);
}
//for int
void save_cat_write_i(hid_t fid, int *data, int ndata, char *fieldname){

	hid_t dspace_id, dset_id;

	//-----
	// Data SPACE
	//-----
	hsize_t dims[1] = {ndata};
	dspace_id = H5Screate_simple(1, dims, NULL);
	
	//-----
	// Set Space
	//-----
	//printf(fieldname);
	dset_id = H5Dcreate(fid, fieldname, H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//-----
	// SAVE
	//-----
	herr_t status;
	status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	//-----
	// Close
	//-----
	H5Dclose(dset_id);
	H5Sclose(dspace_id);
}
//for int64
void save_cat_write_i64(hid_t fid, long long *data, long long ndata, char *fieldname){

	hid_t dspace_id, dset_id;

	//-----
	// Data SPACE
	//-----
	hsize_t dims[1] = {ndata};
	dspace_id = H5Screate_simple(1, dims, NULL);
	
	//-----
	// Set Space
	//-----
	//printf(fieldname);
	dset_id = H5Dcreate(fid, fieldname, H5T_NATIVE_LLONG, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//-----
	// SAVE
	//-----
	herr_t status;
	status = H5Dwrite(dset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	//printf(status);
	//-----
	// Close
	//-----
	H5Dclose(dset_id);
	H5Sclose(dspace_id);
}
//for string
void save_cat_write_s(hid_t fid, char **data, int ndata, char *fieldname){

	hid_t dspace_id, dset_id;
	//-----
	// Data SPACE
	//-----
	size_t max_strlen = 0;
	for (int i=0; i<ndata; i++){
		size_t len = strlen(data[i]);
		if(len > max_strlen) max_strlen = len;
	}

	hid_t str_type = H5Tcopy(H5T_C_S1);

	H5Tset_size(str_type, max_strlen + 1);
	H5Tset_strpad(str_type, H5T_STR_NULLTERM);
	H5Tset_cset(str_type, H5T_CSET_ASCII);

	hsize_t dims[1] = {ndata};
	dspace_id = H5Screate_simple(1, dims, NULL);
	
	//-----
	// Set Space
	//-----
	//printf(fieldname);
	dset_id = H5Dcreate(fid, fieldname, str_type, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//-----
	// SAVE
	//-----
	// String array to contiguous memory
    	char (*flat_data)[max_strlen + 1] = malloc(ndata * (max_strlen + 1));
    	for (int i = 0; i < ndata; i++) {
        	strncpy(flat_data[i], data[i], max_strlen + 1);
    	}

	herr_t status;
	status = H5Dwrite(dset_id, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_data);

	free(flat_data);

	//printf(status);
	//-----
	// Close
	//-----
	H5Dclose(dset_id);
	H5Sclose(dspace_id);
	H5Tclose(str_type);
}
