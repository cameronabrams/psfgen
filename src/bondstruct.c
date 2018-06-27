#include "bondstruct.h"

bondstruct * new_bondstruct ( int * ia, int na ) {
   int i,j;
   bondstruct * bs = (bondstruct*)malloc(sizeof(bondstruct));

   bs->na=na;
   bs->mb=5; // by assumption
   bs->ia=(int*)malloc(na*sizeof(int));
   for (i=0;i<na;i++) bs->ia[i]=ia[i];
   bs->ba=(int**)malloc(na*sizeof(int*));
   for (i=0;i<na;i++) {
     bs->ba[i]=(int*)malloc(4*sizeof(int));
     for (j=0;j<bs->mb;j++) bs->ba[i][j]=-1;
   }
   return bs;
}

void print_bondlist ( bondstruct * bs ) {
   int i,j;
   printf("CFABOND/C) Bondlist:\n");
   for (i=0;i<bs->na;i++) {
      printf("%i : ",bs->ia[i]);
      for (j=0;j<bs->mb&&bs->ba[i][j]!=-1;j++) printf("%i ",bs->ba[i][j]);
      printf("\n");
   }
}

void bondstruct_addbondlist ( bondstruct * bs, int a, int * bl, int nb ) {
  if (bs) {
    int i,ia;
//    printf("addbondlist atomindex %i numbonds %i\n",a,nb);fflush(stdout);
//    printf("searching atomlist\n");fflush(stdout);
    for (ia=0;bs->ia[ia]!=a && ia<bs->na;ia++);
//    printf("result: %i\n",ia);fflush(stdout);
    if (ia==bs->na) {
       printf("ERROR: atom %i is not found in this bondstruct's atomlist\n", a);
       return;
    }
    if (nb>bs->mb) {
       printf("ERROR: atom %i has too many bonds\n",a);
       return;
    } else {
      for (i=0;i<nb;i++) {
         bs->ba[ia][i]=bl[i];
      }
    }
  }
}  

int * bondstruct_getia ( bondstruct * bs ) {
   return bs->ia;
}

int bondstruct_getna ( bondstruct * bs ) {
   return bs->na;
}

int * bondstruct_getrl ( bondstruct * bs, int i, int j ) {
   int * rm = (int*)malloc(bs->na*sizeof(int));
   int i,li,lj, nr=0;
   for (i=0;i<bs->na;i++) rm[i]=-1;

   for (li=0;li<bs->na && bs->ia[li]!=i;li++);
   if (li==bs->na) {
     printf("ERROR: atom %i is not in the bondstruct\n",i);
   }
   for (lj=0;lj<bs->na && bs->ia[lj]!=j;lj++);
   if (lj==bs->na) {
     printf("ERROR: atom %i is not in the bondstruct\n",j);
   }
   rm[nr++]=j;
   for (i=0;i<bs->ba[lj][i]!=-1;i++) rm[nr++]=bs->ba[lj][i];
   grow=1;
   while (grow) {
     grow=0;
     // in progress
   }
   return rm;
}
 
void free_bondstruct ( bondstruct * bs ) {
   int i;
   if (bs->ia) free(bs->ia);
   if (bs->ba) {
      for (i=0;i<bs->na;i++) free(bs->ba[i]);
      free(bs->ba);
   }
}

void free_intarray ( int * a ) {
   free(a);
}
