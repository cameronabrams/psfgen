#ifndef _BONDSTRUCT_H_
#define _BONDSTRUCT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct BONDSTRUCT {
  int na;  // number of atoms
  int mb;  // max neighbors (assumed to be four)
  int * ia; // array of atom indices
  int * rl; // array of atom indices on the "rotation list"
  int nr; // number of atoms on the rotation list
  int sr; // position saver 
  int ** ba; // atom-by-atom array of bond partners, parallel to ia
  int nb; // number of bonds
  int ** b; // list of bonds
} bondstruct;

bondstruct * new_bondstruct ( int * ia, int na );
void free_bondstruct ( bondstruct * bs );
void print_bondlist ( bondstruct * bs );
void bondstuct_makebondlist ( bondstruct * bs );
void bondstruct_addbonds ( bondstruct * bs, int a, int * ba, int nb );
int * bondstruct_getia ( bondstruct * bs );
int bondstruct_getna ( bondstruct * bs );
int bondstruct_getnb ( bondstruct * bs );
int *  bondstruct_getrl ( bondstruct * bs, int a, int b );
int bondstruct_arebonded ( bondstruct * bs, int a, int b );
void free_intarray ( int * a );
int * bondstruct_getbondpointer ( bondstruct * bs, int i );

#endif
