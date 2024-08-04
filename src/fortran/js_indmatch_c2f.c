#include <stdio.h>

void js_indmatch(int argc, void *argv[])
{
  extern void js_indmatch_();   /* FORTRAN routine */
  int *larr;
  double *darr;
  long long *x, *y, *x_match, *y_match;

  larr  = (int *) argv[0];
  darr  = (double *) argv[1];
  x		= (long long *) argv[2];  
  y		= (long long *) argv[3];
  x_match	= (long long *) argv[4];
  y_match	= (long long *) argv[5];  

  js_indmatch_(larr, darr, x, y, x_match, y_match);   /* Compute sum */
}
