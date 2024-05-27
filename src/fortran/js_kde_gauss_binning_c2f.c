#include <stdio.h>

void js_kde_gauss_binning(int argc, void *argv[])
{
  extern void js_kde_gauss_binning_();   /* FORTRAN routine */
  double *xx, *yy, *ww;
  double *grid;
  double *xr, *yr;
  int *larr;

  xx	= (double *) argv[0];  
  yy	= (double *) argv[1];
  ww	= (double *) argv[2];
  xr	= (double *) argv[3];
  yr	= (double *) argv[4];
  grid	= (double *) argv[5];
  larr	= (int *) argv[6];

  js_kde_gauss_binning_(xx, yy, ww, xr, yr, grid, larr);   /* Compute sum */
}
