#include <stdio.h>

void js_kde_2dhisto(int argc, void *argv[])
{
  extern void js_kde_2dhisto_();   /* FORTRAN routine */
  double *xx, *yy, *ww;
  float *grid, *grid2;
  double *xr, *yr;
  int *n_pix;
  int *larr;
  double *darr;

  xx	= (double *) argv[0];  
  yy	= (double *) argv[1];
  ww	= (double *) argv[2];
  grid	= (float *) argv[3];
  grid2 = (float *) argv[4];
  xr	= (double *) argv[5];
  yr	= (double *) argv[6];
  larr	= (int *) argv[7];

  js_kde_2dhisto_(xx, yy, ww, grid, grid2, xr, yr, larr);   /* Compute sum */
}
