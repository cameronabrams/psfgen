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
  int nb; // number of rotatable bonds
  int ** b; // array of pairs of atom indices of rotatable bonds
  int ** bra; // array of right-side atom indicies for each rotatable bond
  int * bran; // array of counts of right-side atoms for each rotatable bond
  int * isactive; // array of is-active flags; 0 means bond is not active for rotation
} bondstruct;

bondstruct * new_bondstruct ( int * ia, int na );
void free_bondstruct ( bondstruct * bs );
void print_bondlist ( bondstruct * bs );
void bondstruct_makerotatablebondlist ( bondstruct * bs, int * rot_i, int ni, int * rot_j, int nj );
void bondstruct_addbonds ( bondstruct * bs, int a, int * ba, int nb );
int * bondstruct_getia ( bondstruct * bs );
int bondstruct_getna ( bondstruct * bs );
int bondstruct_getnb ( bondstruct * bs );
int bondstruct_arebonded ( bondstruct * bs, int a, int b );
void free_intarray ( int * a );
int * bondstruct_getbondpointer ( bondstruct * bs, int i );
int * bondstruct_getrightside_pointer ( bondstruct * bs, int b );
int bondstruct_getrightside_count ( bondstruct * bs, int b );
int bondstruct_getbondindex ( bondstruct * bs, int i, int j );
double my_roughenergy ( int * i1, double * x1, double * y1, double * z1, int n1, int * i2, 
double * x2, double * y2, double * z2, int n2, double cut, double sigma, double epsilon, 
bondstruct * bs );
void bondstruct_deactivate_by_fixed ( bondstruct * bs, int fa );
int bondstruct_isactive ( bondstruct * bs, int b ); 
#endif
