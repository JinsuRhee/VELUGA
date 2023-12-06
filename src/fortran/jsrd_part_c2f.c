#include <stdio.h>

typedef struct {
	unsigned short slen;
	short stype;
	char *s;
} STRING;

#define STR_LEN(__str)	((long)(__str)->slen)

void jsrd_part(int argc, void *argv[])
{
	extern void jsrd_part_();
	int *larr;
	double *darr;
	double *xp, *vp, *mp, *ap, *zp;
	int *fam, *tag, *domain;
	int *part_ind, *domlist;
	long long *idvar;
	STRING *file;

	larr	= (int *) argv[0];
	darr	= (double *) argv[1];
	file	= (STRING *) argv[2];
	part_ind= (int *) argv[3];
	xp	= (double *) argv[4];
	vp	= (double *) argv[5];
	mp	= (double *) argv[6];
	ap	= (double *) argv[7];
	zp	= (double *) argv[8];
	fam	= (int *) argv[9];
	tag	= (int *) argv[10];
	domain	= (int *) argv[11];
	idvar	= (long long *) argv[12];
	domlist	= (int *) argv[13];

	jsrd_part_(larr, darr, file->s, part_ind, xp, vp, mp, ap, zp, fam, tag, domain, idvar, domlist);
}
	
