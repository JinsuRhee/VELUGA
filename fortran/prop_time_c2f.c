#include <stdio.h>

void prop_time(int argc, void *argv[])
{
  extern void prop_time_();   /* FORTRAN routine */
  double *p_age, *p_sf, *p_gyr;
  double *tmp_c, *tmp_s, *tmp_r, *tmp_g;
  int *larr;
  double *darr;

  p_age		= (double *) argv[0];
  p_sf		= (double *) argv[1];
  p_gyr		= (double *) argv[2];
  tmp_c		= (double *) argv[3];
  tmp_s		= (double *) argv[4];
  tmp_r		= (double *) argv[5];
  tmp_g		= (double *) argv[6]; 
  larr		= (int *) argv[7];
  darr		= (double *) argv[8];

  prop_time_(p_age, p_sf, p_gyr, tmp_c, tmp_s, tmp_r, tmp_g, larr, darr);   /* Compute sum */
}
