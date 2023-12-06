#include <stdio.h>
#include <sys/resource.h>

typedef struct {
   unsigned short slen;         /* length of the string         */
   short stype;                 /* Type of string               */
   char *s;                     /* Pointer to chararcter array  */
} STRING;

#define STR_LEN(__str)    ((long)(__str)->slen)

void get_contam(int argc, void *argv[])
{
  extern void get_contam_();   /* FORTRAN routine */
  double *xc, *yc, *zc, *rr;
  int *larr;
  double *darr;
  double *dm_xp, *dm_mm;
  double *conf_n, *conf_m, *conf_r;

  STRING *dir_raw;

  larr		= (int *) argv[0];
  darr		= (double *) argv[1];
  dir_raw	= (STRING *) argv[2];
  xc		= (double *) argv[3];
  yc		= (double *) argv[4];
  zc		= (double *) argv[5];
  rr		= (double *) argv[6];
  dm_xp		= (double *) argv[7];
  dm_mm		= (double *) argv[8];
  conf_n	= (double *) argv[9];
  conf_m	= (double *) argv[10];
  conf_r	= (double *) argv[11];


  get_contam_(larr, darr, dir_raw->s, xc, yc, zc, rr, dm_xp, dm_mm, conf_n, conf_m, conf_r);   /* Compute sum */}
