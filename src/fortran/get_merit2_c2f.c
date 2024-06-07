#include <stdio.h>

void get_merit2(int argc, void *argv[])
{
  extern void get_merit2_();   /* FORTRAN routine */
  int *larr;
  double *darr;
  long long *pid_g, *pid_s;
  long *gid_g, *gid_s, *npart_s, *npart_g;
  long *hash, *hash_next, *m_id;
  double *merit, *m_merit;

  larr		= (int *) argv[0];
  darr		= (double *) argv[1];
  pid_g		= (long long *) argv[2];
  gid_g		= (long *) argv[3];
  pid_s		= (long long *) argv[4];
  gid_s		= (long *) argv[5];
  hash		= (long *) argv[6];
  hash_next	= (long *) argv[7];
  npart_g	= (long *) argv[8];
  npart_s	= (long *) argv[9];
  merit		= (double *) argv[10];
  m_id		= (long *) argv[11];
  m_merit	= (double *) argv[12];

  get_merit2_(larr, darr, pid_g, gid_g, pid_s, gid_s, hash, hash_next, npart_g, npart_s, merit, m_id, m_merit);   /* Compute sum */
}
