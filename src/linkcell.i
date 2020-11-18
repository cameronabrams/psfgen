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
#include "linkcell.h"
%}
extern linkcell * linkcell_new ( double * x, double * y, double * z, int n, double cut );
extern linkcell * my_roughenergy_setup ( double * x2, double * y2, double * z2, int n2, double cut );
