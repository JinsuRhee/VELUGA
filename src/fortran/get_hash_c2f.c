#include <stdio.h>

void get_hash(int argc, void *argv[])
{
  extern void get_hash_();   /* FORTRAN routine */
  int *larr;
  double *darr;
  long long *pid;
  long *indarr, *hash, *hash_next;

  larr		= (int *) argv[0];
  darr		= (double *) argv[1];
  pid		= (long long *) argv[2];
  indarr	= (long *) argv[3];
  hash		= (long *) argv[4];
  hash_next	= (long *) argv[5];

  get_hash_(larr, darr, pid, indarr, hash, hash_next);   /* Compute sum */
}
