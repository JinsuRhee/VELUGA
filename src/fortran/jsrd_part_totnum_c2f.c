#include <stdio.h>

typedef struct {
	unsigned short slen;
	short stype;
	char *s;
} STRING;

#define STR_LEN(__str)	((long)(__str)->slen)

void jsrd_part_totnum(int argc, void *argv[])
{
	extern void jsrd_part_totnum_();
	int *larr, *npart_tot, *part_ind, *domlist;
	double *darr;
	STRING *file;

	larr	= (int *) argv[0];
	darr	= (double *) argv[1];
	file	= (STRING *) argv[2];
	npart_tot=(int *) argv[3];
	part_ind= (int *) argv[4];
	domlist	= (int *) argv[5];

	jsrd_part_totnum_(larr, darr, file->s, npart_tot, part_ind, domlist);
}
	
