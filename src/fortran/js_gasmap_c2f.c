#include <stdio.h>

void js_gasmap(int argc, void *argv[])
{
  extern void js_gasmap_();   /* FORTRAN routine */
  int *larr;
  double *darr;
  double *xx, *yy, *zz, *bw, *xr, *yr, *map;

  larr	= (int *) 	argv[0];
  darr	= (double *)	argv[1];
  xx	= (double *)	argv[2];
  yy	= (double *)	argv[3];
  zz	= (double *)	argv[4];
  bw	= (double *)	argv[5];
  xr	= (double *)	argv[6];
  yr	= (double *)	argv[7];
  map	= (double *)	argv[8];

  js_gasmap_(larr, darr, xx, yy, zz, bw, xr, yr, map);   /* Compute sum */
}
