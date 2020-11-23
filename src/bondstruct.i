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
#include "linkcell.h"
%}
extern bondstruct * new_bondstruct ( int * ia, int na, int nb );
extern void free_bondstruct ( bondstruct * bs );
extern void bondstruct_print ( bondstruct * bs );
extern int bondstruct_getlocalatomindex ( bondstruct * bs, int a );
extern int bondstruct_getbondindex ( bondstruct * bs, int i, int j );
extern int * bondstruct_getrightside_pointer ( bondstruct * bs, int b );
extern int bondstruct_getrightside_count ( bondstruct * bs, int b );
extern void bondstruct_addbond( bondstruct * bs, int i, int j);
extern void bondstruct_importbonds ( bondstruct * bs, int a, int * bl, int nb );
extern void bondstruct_setbond_rotatable ( bondstruct * bs, int a, int b, int flag );
extern void bondstruct_setbond_active ( bondstruct * bs, int a, int b, int flag );
extern void bondstruct_maprotatables ( bondstruct * bs );
extern int bondstruct_r2b ( bondstruct * bs, int r );
extern int bondstruct_getnb ( bondstruct * bs );
extern int bondstruct_getnrb ( bondstruct * bs );
extern int * bondstruct_getbondpointer ( bondstruct * bs, int i );
extern int * bondstruct_getia ( bondstruct * bs );
extern int bondstruct_getna ( bondstruct * bs );
extern void bondstruct_makerightsides ( bondstruct * bs );
extern int bondstruct_isactive ( bondstruct * bs, int b );
extern void bondstruct_deactivate_by_fixed ( bondstruct * bs, int fa );
extern int bondstruct_arebonded ( bondstruct * bs, int a, int b );
extern linkcell * my_roughenergy_setup ( double * x2, double * y2, double * z2, int n2, double cut );
extern double my_roughenergy ( int * i1, double * x1, double * y1, double * z1, int n1, int * i2, 
                        int n2, double cutoff,
                        double sigma, double epsilon, double shift,
                        bondstruct * bs, linkcell * ls );
extern void my_roughenergy_cleanup ( linkcell * ls );
