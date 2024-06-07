#include <stdio.h>

void get_merit(int argc, void *argv[])
{
  extern void get_merit_();   /* FORTRAN routine */
  int *larr;
  double *darr;
  double *w0, *w1;
  long long *pid, *pid2;

  larr		= (int *) argv[0];
  darr		= (double *) argv[1];
  pid		= (long long *) argv[2];
  pid2		= (long long *) argv[3];
  w0		= (double *) argv[4];
  w1		= (double *) argv[5];

  get_merit_(larr, darr, pid, pid2, w0, w1);   /* Compute sum */
}
