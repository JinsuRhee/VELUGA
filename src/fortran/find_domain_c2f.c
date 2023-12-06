#include <stdio.h>

void find_domain(int argc, void *argv[])
{
  extern void find_domain_();   /* FORTRAN routine */
  double *xc, *yc, *zc, *rr, *hindex;
  long *levmax, *dom_list;
  long *larr;
  double *darr;

  xc		= (double *) argv[0];
  yc		= (double *) argv[1];
  zc		= (double *) argv[2];
  rr		= (double *) argv[3];
  hindex	= (double *) argv[4];
  levmax	= (long *) argv[5];
  dom_list	= (long *) argv[6];
  larr		= (long *) argv[7];
  darr		= (double *) argv[8];

  find_domain_(xc, yc, zc, rr, hindex, levmax, dom_list, larr, darr);   /* Compute sum */
}
