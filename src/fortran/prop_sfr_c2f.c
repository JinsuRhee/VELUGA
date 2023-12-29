#include <stdio.h>

void prop_sfr(int argc, void *argv[])
{
  extern void prop_sfr_();   /* FORTRAN routine */
  double *xc, *yc, *zc;
  double *xx, *yy, *zz, *age, *mass;
  int *larr;
  double *darr;

  larr		= (int *) argv[0];
  darr		= (double *) argv[1];
  xx		= (double *) argv[2];
  yy	 	= (double *) argv[3];
  zz		= (double *) argv[4];
  age		= (double *) argv[5];
  mass		= (double *) argv[6];
  xc		= (double *) argv[7];
  yc		= (double *) argv[8];
  zc		= (double *) argv[9];


  prop_sfr_(larr, darr, xx, yy, zz, age, mass, xc, yc, zc);   /* Compute sum */
}
