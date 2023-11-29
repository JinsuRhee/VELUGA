#include <stdio.h>

typedef struct {
   unsigned short slen;         /* length of the string         */
   short stype;                 /* Type of string               */
   char *s;                     /* Pointer to chararcter array  */
} STRING;

#define STR_LEN(__str)    ((long)(__str)->slen)

void get_ptcl(int argc, void *argv[])
{
  extern void get_ptcl_();   /* FORTRAN routine */
  long long *id;
  double *ptcl;
  long *dom_list;

  STRING *dir_raw; 
  int *larr;
  double *darr;

  larr		= (int *) argv[0];
  darr		= (double *) argv[1];

  dir_raw	= (STRING *) argv[2];
  id		= (long long *) argv[3];
  ptcl		= (double *) argv[4];
  dom_list	= (long *) argv[5];

  //printf("%s", dir_raw->s);

  get_ptcl_(larr, darr, dir_raw->s, id, ptcl, dom_list);   /* Compute sum */
}
