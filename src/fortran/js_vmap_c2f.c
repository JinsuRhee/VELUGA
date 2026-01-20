#include <stdio.h>

void js_vmap(int argc, void *argv[])
{
  extern void js_vmap_();   /* FORTRAN routine */
  int *larr;
  double *darr;
  int *ix, *iy;
  double *vx, *vy, *mass, *vmap;

  larr	= (int *) 	argv[0];
  darr	= (double *)	argv[1];
  ix	= (int *)	argv[2];
  iy	= (int *)	argv[3];
  vx	= (double *)	argv[4];
  vy	= (double *)	argv[5];
  mass	= (double *)	argv[6];
  vmap	= (double *)	argv[7];


  js_vmap_(larr, darr, ix, iy, vx, vy, mass, vmap);   /* Compute sum */
}
