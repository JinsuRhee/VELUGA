#include <stdio.h>

typedef struct {
	unsigned short slen;
	short stype;
	char *s;
} STRING;

#define STR_LEN(__str)	((long)(__str)->slen)

void jsamr2cell(int argc, void *argv[])
{
	extern void jsamr2cell_();
	int *larr, *mg_ind, *mesh_lv, *domlist;
	double *darr, *mesh_xg, *mesh_hd, *mesh_dx;
	STRING *file_a, *file_h, *file_i;

	larr	= (int *) argv[0];
	darr	= (double *) argv[1];
	file_a	= (STRING *) argv[2];
	file_h	= (STRING *) argv[3];
	file_i	= (STRING *) argv[4];
	mg_ind	= (int *) argv[5];
	mesh_xg	= (double *) argv[6];
	mesh_dx	= (double *) argv[7];
	mesh_hd	= (double *) argv[8];
	mesh_lv	= (int *) argv[9];
	domlist = (int *) argv[10];

	jsamr2cell_(larr, darr, file_a->s, file_h->s, file_i->s, mg_ind, mesh_xg, mesh_dx, mesh_hd, mesh_lv, domlist);
}
	
