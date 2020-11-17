%module cfa_bondstruct
%include carrays.i
%array_functions(double, array);
%array_functions(int, arrayint);
%inline %{
double get_double(double *a, int index) {
  return a[index];
}
%}
%{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bondstruct.h"
%}
extern bondstruct * new_bondstruct ( int * ia, int na );
extern void print_bondlist ( bondstruct * bs );
extern void bondstruct_addbonds ( bondstruct * bs, int a, int * ba, int nb );
extern void bondstruct_makerotatablebondlist ( bondstruct * bs, int * rot_i, int ni, int * rot_j, int nj );
extern void free_bondstruct ( bondstruct * bs );
extern int bondstruct_getna ( bondstruct * bs );
extern int bondstruct_getnb ( bondstruct * bs );
extern int bondstruct_arebonded ( bondstruct * bs, int a, int b );
extern int * bondstruct_getbondpointer ( bondstruct * bs, int i );
extern int * bondstruct_getrightside_pointer ( bondstruct * bs, int b );
extern int bondstruct_getrightside_count ( bondstruct * bs, int b );
extern int bondstruct_getbondindex ( bondstruct * bs, int i, int j );

extern void free_intarray ( int * a );
extern double my_roughenergy ( int * i1, double * x1, double * y1, double * z1, int n1, int * i2, 
double * x2, double * y2, double * z2, int n2, double cut, double sigma, double epsilon );
extern void bondstruct_deactivate_by_fixed ( bondstruct * bs, int fa );
extern int bondstruct_isactive ( bondstruct * bs, int b ); 

