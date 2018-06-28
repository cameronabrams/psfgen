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
extern void bondstruct_makebondlist ( bondstruct * bs );
extern void free_bondstruct ( bondstruct * bs );
extern int * bondstruct_getrl ( bondstruct * bs, int a, int b );
extern int bondstruct_getna ( bondstruct * bs );
extern int bondstruct_getnb ( bondstruct * bs );
extern int bondstruct_arebonded ( bondstruct * bs, int a, int b );
extern int * bondstruct_getbondpointer ( bondstruct * bs, int i );
extern void free_intarray ( int * a );
