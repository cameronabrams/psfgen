#include "bondstruct.h"

bondstruct * new_bondstruct ( int * ia, int na ) {
   int i,j;
   bondstruct * bs = (bondstruct*)malloc(sizeof(bondstruct));
   bs->nr=0;
   bs->sr=0;
   bs->na=na;
   bs->mb=5; // by assumption
   bs->ia=(int*)malloc(na*sizeof(int));
   for (i=0;i<na;i++) bs->ia[i]=ia[i];
   bs->rl=(int*)malloc(na*sizeof(int));
   for (i=0;i<na;i++) bs->rl[i]=0;
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


int bondstruct_getlocalindex ( bondstruct * bs, int a ) {
   int la;
   printf("%i\n",a);fflush(stdout);
   for (la=0;la<bs->na && bs->ia[la]!=a;la++);
   if (la==bs->na) {
     printf("ERROR: atom %i is not in the bondstruct\n",a);
     return -1;
   }
   return la;
}

void bondstruct_resetrotationlist ( bondstruct * bs ) {
   int i;   
   bs->nr=bs->sr=0;
   for (i=0;i<bs->na;i++) bs->rl[i]=-1;
}

// set the rotation list based on the bond arrays and the identities
// of the two bonded atoms; "a" is the upstream atom and "b" the
// downstream.  Atom "b" and the rest of the molecule accessible
// "b" via bond traversal excluding the bond to atom "a" is rotatable 
int * bondstruct_getrl ( bondstruct * bs, int a, int b ) {
   int i,j,k,l,m,la,lb,nr=0,grow,lnr,ca,lk;
   bondstruct_resetrotationlist(bs);
   printf("%i %i\n",a,b); fflush(stdout);
   lb=bondstruct_getlocalindex(bs,b);
   // put the downstream atom at the head of the rotation list
   bs->rl[bs->nr++]=b;
   bs->sr++;
   // put all atoms it is bonded to on the rotation list EXCEPT the upstream!
   for (i=0;i<bs->mb && bs->ba[lb][i]!=-1;i++) {
     if (bs->ba[lb][i] != a) {
       bs->rl[bs->nr++]=bs->ba[lb][i];
     }
   }
   // grow
   grow=1;
   while (grow) {
     grow=0;
     // between sr and nr-1, there are atoms whose neighbors should be added to rotlist
     lnr=bs->nr;
     for (k=bs->sr;k<lnr;k++) {
        lk=bondstruct_getlocalindex(bs,bs->rl[k]);
        for (l=0;l<bs->mb && bs->ba[lk][l]!=-1;l++) {
          if (bs->ba[lk][l] != a) {
            // only add this to the rotlist if it is not already on the rotlist
            ca=bs->ba[lk][l];
            for (m=0;m<bs->nr&&bs->rl[m]!=ca;m++);
            if (ca!=bs->rl[m]) {
              bs->rl[bs->nr++]=bs->ba[lk][l];
              grow=1;
            }
          }
        }
     }
     if (grow) bs->sr=lnr;
   }
   return bs->rl;
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
