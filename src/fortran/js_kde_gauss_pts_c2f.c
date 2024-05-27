#include <stdio.h>

void js_kde_gauss_pts(int argc, void *argv[])
{
  extern void js_kde_gauss_pts_();   /* FORTRAN routine */
  double *xx, *yy, *ww, *ptcl;
  double *grid;
  double *xr, *yr;
  int *larr;

  xx	= (double *) argv[0];  
  yy	= (double *) argv[1];
  ww	= (double *) argv[2];
  xr	= (double *) argv[3];
  yr	= (double *) argv[4];
  grid	= (double *) argv[5];
  ptcl	= (double *) argv[6];
  larr	= (int *) argv[7];

  js_kde_gauss_pts_(xx, yy, ww, xr, yr, grid, ptcl, larr);   /* Compute sum */
}
