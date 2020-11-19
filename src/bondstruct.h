#ifndef _BONDSTRUCT_H_
#define _BONDSTRUCT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linkcell.h"

// rotatable bonds: bonds are directional!  the i->j bond is "different"
// from the j->i bond because even though the same two atoms
// are involved, rotating the i->j bond moves atoms in the part
// bound to j; i.e., j is the "right-atom" of the i->j bond
// while rotating the j->i bond would move all atoms in the
// part bound to i; i.e., i is the "right-atom" of the j->i bond
//
// [$sel getbonds] in vmd returns a lists of lists.  Each element
// is the list of atom indicies of atoms bound to the atom that owns
// that element.  The order of the elements is the same as the order
// of atom indicies in [$sel get index].  The element lists contain
// all bonded partners; if i and j are bonded, then j appears on i's 
// list and i appears on j's list.
//
// So *every* bond in the bondstruct that is rotatable has a "right-side"
// which is all atoms in the sel that would move if the bond is rotated
typedef struct BONDSTRUCT {
  int na;  // number of atoms
  int mb;  // max neighbors (assumed to be four)
  int * ia; // array of atom indices dim [na]
  int ** ba; // atom-by-atom array of bond partners dim [na][mb]
  int nb; // number of bonds in bondlist (all bonds)
  int * isrotatable; // array of is-rotatable flags; 1 means bond is rotatable dim [nb]
  // the four arrays below only apply for rotatable bonds, but every bond gets an entry
  int * isactive; // array of is-active flags; 0 means bond is not active for rotation dim [nb]
  int * bran; // array of counts of right-side atoms for each rotatable bond dim [nb]
  int bctr;
  int ** b; // array of pairs of atom indices of rotatable bonds dim [nb][2]
  int ** bra; // array of right-side atom indicies for each rotatable bond dim [nb][bran[i]]

  int nrb; // number of rotatable bonds, equals number of 1's in isrotatable[]
  int * r2b; // maps index [0:nrb-1] to index [0:nb-1], dim[nrb]

  //int nr; // number of atoms on the rotation list
  //int * rl; // array of atom indices on the "rotation list" dim [nr]
  //int sr; // position saver 
} bondstruct;


bondstruct * new_bondstruct ( int * ia, int na, int nb );
void free_bondstruct ( bondstruct * bs );
void bondstruct_print ( bondstruct * bs );
int bondstruct_getlocalatomindex ( bondstruct * bs, int a );
int bondstruct_getbondindex ( bondstruct * bs, int i, int j );
int * bondstruct_getrightside_pointer ( bondstruct * bs, int b );
int bondstruct_getrightside_count ( bondstruct * bs, int b );
void bondstruct_addbond( bondstruct * bs, int i, int j);
void bondstruct_importbonds ( bondstruct * bs, int a, int * bl, int nb );
void bondstruct_setbond_rotatable ( bondstruct * bs, int a, int b, int flag );
void bondstruct_setbond_active ( bondstruct * bs, int a, int b, int flag );
void bondstruct_maprotatables ( bondstruct * bs );
int bondstruct_r2b ( bondstruct * bs, int r );
int bondstruct_getnb ( bondstruct * bs );
int bondstruct_getnrb ( bondstruct * bs );
int * bondstruct_getbondpointer ( bondstruct * bs, int i );
int * bondstruct_getia ( bondstruct * bs );
int bondstruct_getna ( bondstruct * bs );
void bondstruct_makerightsides ( bondstruct * bs );
int bondstruct_isactive ( bondstruct * bs, int b );
void bondstruct_deactivate_by_fixed ( bondstruct * bs, int fa );
int bondstruct_arebonded ( bondstruct * bs, int a, int b );

double my_roughenergy ( int * i1, double * x1, double * y1, double * z1, int n1, 
                        int * i2, int n2, double cutoff,
                        double sigma, double epsilon, double shift,
                        bondstruct * bs, linkcell * ls );
void my_roughenergy_cleanup ( linkcell * ls );
#endif
