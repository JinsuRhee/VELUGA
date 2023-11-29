#include <stdio.h>

typedef struct {
	unsigned short slen;
	short stype;
	char *s;
} STRING;

#define STR_LEN(__str)	((long)(__str)->slen)

void jsamr2cell_totnum(int argc, void *argv[])
{
	extern void jsamr2cell_totnum_();
	int *larr, *ntot, *nvarh, *mg_ind, *domlist;
	double *darr;
	STRING *file_a, *file_h;

	larr	= (int *) argv[0];
	darr	= (double *) argv[1];
	file_a	= (STRING *) argv[2];
	file_h	= (STRING *) argv[3];
	ntot	= (int *) argv[4];
	nvarh	= (int *) argv[5];
	mg_ind	= (int *) argv[6];
	domlist	= (int *) argv[7];

	jsamr2cell_totnum_(larr, darr, file_a->s, file_h->s, ntot, nvarh, mg_ind, domlist);
}
	
