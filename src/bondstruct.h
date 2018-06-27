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
  int ** ba; // array of bond partners
} bondstruct;

bondstruct * new_bondstruct ( int * ia, int na );
void free_bondstruct ( bondstruct * bs );
void print_bondlist ( bondstruct * bs );
void bondstruct_addbondlist ( bondstruct * bs, int a, int * ba, int nb );
int * bondstruct_getia ( bondstruct * bs );
int bondstruct_getna ( bondstruct * bs );
int *  bondstruct_getrl ( bondstruct * bs, int a, int b );
void free_intarray ( int * a );

#endif
