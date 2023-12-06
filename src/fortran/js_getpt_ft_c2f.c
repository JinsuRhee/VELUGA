#include <stdio.h>
#include <sys/resource.h>

void js_getpt_ft(int argc, void *argv[])
{
  extern void js_getpt_ft_();   /* FORTRAN routine */
  int *larr;
  double *darr;

  double *p_xp, *p_mp, *p_pot, *p_force;

  larr		= (int *) argv[0]; 
  darr		= (double *) argv[1];
  p_xp		= (double *) argv[2];
  p_mp		= (double *) argv[3];
  p_pot		= (double *) argv[4];
  p_force	= (double *) argv[5];

  js_getpt_ft_(larr, darr, p_xp, p_mp, p_pot, p_force);   /* Compute sum */}
