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
   bs->b=NULL;
   bs->nb=0;
   bs->bra=NULL;
   bs->bran=NULL;
   return bs;
}

void free_bondstruct ( bondstruct * bs ) {
   int i;
   if (bs->ia) free(bs->ia);
   if (bs->rl) free(bs->rl);
   if (bs->ba) {
      for (i=0;i<bs->na;i++) free(bs->ba[i]);
      free(bs->ba);
   }
   if (bs->b) free(bs->b);
   if (bs->bra) free(bs->bra);
   if (bs->bran) free(bs->bran);
   free(bs);
}

int bondstruct_getlocalindex ( bondstruct * bs, int a ) {
   int la;
   for (la=0;la<bs->na && bs->ia[la]!=a;la++);
   if (la==bs->na) {
     printf("ERROR: atom %i is not in the bondstruct\n",a);
     return -1;
   }
   return la;
}

void print_bondlist ( bondstruct * bs ) {
   int i,j;
   printf("CFABOND/C) Array of rotatable bonds:\n");
   if (bs->b && bs->nb) {
      for (i=0;i<bs->nb;i++) {
        printf("%c%i %i : %i : ",(bs->isactive[i]?'+':'-'),bs->b[i][0],bs->b[i][1],bs->bran[i]);
        for (j=0;j<bs->bran[i];j++) printf("%i ",bs->bra[i][j]);
        printf("\n");
      }
   }
   printf("CFABOND/C) Atom bondlists:\n");
   for (i=0;i<bs->na;i++) {
      printf("%i : ",bs->ia[i]);
      for (j=0;j<bs->mb&&bs->ba[i][j]!=-1;j++) printf("%i ",bs->ba[i][j]);
      printf("\n");
   }
}

int bondstruct_getbondindex ( bondstruct * bs, int i, int j ) {
   int k=0;
   int m=-1;
   for (k=0;k<bs->nb;k++) {
     if (bs->b[k][0]==i&&bs->b[k][1]==j) m=k;
   }
   return m;
}

int * bondstruct_getrightside_pointer ( bondstruct * bs, int b ) {
   return bs->bra[b];
}

int bondstruct_getrightside_count ( bondstruct * bs, int b ) {
   return bs->bran[b];
}

// a is an atom index and bl[] is an array of atom indices of its bonding partners
// as determined by [$sel getbonds].  nb is the number of such partners.
void bondstruct_addbonds ( bondstruct * bs, int a, int * bl, int nb ) {
  if (bs) {
    int i,ia;
    ia=bondstruct_getlocalindex(bs,a);
    if (nb>bs->mb) {
       printf("ERROR: atom %i has too many bonds\n",a);
       return;
    } else {
      for (i=0;i<nb;i++) {
         if (bl[i]==0) {
            printf("ERROR: position %i of atom %i neighbor list is zero!\n",i,a);
         }
         bs->ba[ia][i]=bl[i];
      }
    }
  }
}  

int bondstruct_getnb ( bondstruct * bs ) {
   return bs->nb;
}

int * bondstruct_getbondpointer ( bondstruct * bs, int i ) {
   return bs->b[i];
}

int bond_rotatable ( int a, int b, int * ri, int ni, int * rj, int nj ) {
    int i,j;
    // is a on ri?
    for (i=0;i<ni&&ri[i]!=a;i++);
    // is b on rj?
    for (j=0;j<nj&&rj[j]!=b;j++);
    if (i<ni&&j<nj) return 1;
    return 0;
} 

int * bondstruct_getia ( bondstruct * bs ) {
   return bs->ia;
}

int bondstruct_getna ( bondstruct * bs ) {
   return bs->na;
}


// After making the rotatable bond list, this function 
// is called to make, for _each_ rotatable bond, the
// array of atom indices that are on the "right" side 
// of the bond. These are atoms that move when the 
// bond is later rotated, and they are a static 
// property of the bondstruct that need be determined 
// only once
void bondstruct_makerightsides ( bondstruct * bs ) {
   int i,j,k,l,m,a,b,la,lb,li,nr=0,grow,lnr,ca,lk;
   
   // allocate the array of arrays
   bs->bra=(int**)malloc(bs->nb*sizeof(int*));
   for (i=0;i<bs->nb;i++) {
     bs->bra[i]=(int*)malloc(bs->na*sizeof(int));
     for (j=0;j<bs->na;j++) bs->bra[i][j]=-1;
   }
   // allocate array of right-side counts
   bs->bran=(int*)malloc(bs->nb*sizeof(int));
   // for each bond, identify and save all right-side atoms
   //printf("#### in makerightsides processing %i bonds\n",bs->nb);fflush(stdout);
   for (k=0;k<bs->nb;k++) {
     a=bs->b[k][0];
     b=bs->b[k][1];
     //printf("#### in makerightsides at bond %i : %i %i\n",k,a,b);fflush(stdout);
     la=bondstruct_getlocalindex(bs,a);
     lb=bondstruct_getlocalindex(bs,b);
     //printf("#### local indices %i %i\n",la,lb);fflush(stdout);
     // put all atoms b is bonded to on the right-side list _except_ atom a
     bs->bran[k]=0;
     for (i=0;i<bs->mb && bs->ba[lb][i]!=-1;i++) {
       if (bs->ba[lb][i] != a) {
         //printf("#### added %i-neighbor %i to bond-%i rightside list\n",b,bs->ba[lb][i],k);fflush(stdout);
         bs->bra[k][bs->bran[k]++]=bs->ba[lb][i];
       }
     }
     // enter the grow-loop; for each element on bs->bra[k][], add _its_ bond partners to this
     // right-side array, provided they are not already in it
     //printf("#### growing rightside of bond %i\n",k);fflush(stdout);
     grow=1;
     while (grow) {
       grow=0;
       // save the current right-side array size
       // lnr=bs->bran[k];
       // for each atom on this right-side array
       for (i=0;i<bs->bran[k];i++) {
          // get the local index
          li=bondstruct_getlocalindex(bs,bs->bra[k][i]);
          // for each neighbor atom bonded to this atom
          for (j=0;j<bs->mb && bs->ba[li][j]!=-1;j++) {
            // as long as this neighbor is not atom a (which it should never be)
            if (bs->ba[li][j] != a) {
              // only add this to the right-side array if it is not already in it
              ca=bs->ba[li][j];
              for (m=0;m<bs->bran[k]&&bs->bra[k][m]!=ca;m++);
              if (ca!=bs->bra[k][m]) {
                bs->bra[k][bs->bran[k]++]=ca;
                grow=1;
              }
            }
          }
        }
      }
   }
}

void bondstruct_makerotatablebondlist ( bondstruct * bs, int * rot_i, int ni, int * rot_j, int nj ) {
   int i,j,k;
   bs->nb=0;
   for (i=0;i<bs->na;i++) {
     for (j=0;j<bs->mb&&bs->ba[i][j]!=-1;j++) 
       if (bond_rotatable(bs->ia[i],bs->ba[i][j],rot_i,ni,rot_j,nj)) bs->nb++;
   }
   bs->b=(int**)malloc(bs->nb*sizeof(int*));
   bs->isactive=(int*)malloc(bs->nb*sizeof(int));
   k=0;
   for (i=0;i<bs->na;i++) {
     for (j=0;j<bs->mb&&bs->ba[i][j]!=-1;j++) {
       if (bond_rotatable(bs->ia[i],bs->ba[i][j],rot_i,ni,rot_j,nj)) {
          bs->b[k]=(int*)malloc(2*sizeof(int));
          bs->b[k][0]=bs->ia[i];
          bs->b[k][1]=bs->ba[i][j];
          bs->isactive[k]=1;
          // generate array of atoms on the "right" of this bond
          k++;
       }
     }
   }  
//   printf("#### calling bondstruct_makerightsides\n");fflush(stdout);
   bondstruct_makerightsides(bs);
}

int bondstruct_isactive ( bondstruct * bs, int b ) {
   return bs->isactive[b] && bs->bran[b]>0;
}

void bondstruct_deactivate_by_fixed ( bondstruct * bs, int fa ) {
   int i,j;
   for (i=0;i<bs->nb;i++) {
      for (j=0;j<bs->bran[i]&&bs->bra[i][j]!=fa;j++);
      if (j<bs->bran[i]) bs->isactive[i]=0;
   }
}

int bondstruct_arebonded ( bondstruct * bs, int a, int b ) {
  int rm=0;

  int * ba=bs->ba[bondstruct_getlocalindex(bs,a)];
  int i;
  for (i=0;i<bs->mb&&ba[i]!=-1&&ba[i]!=b;i++);
  if (i==bs->mb||ba[i]==-1) return 0;
  else return rm;
}
 
void free_intarray ( int * a ) {
   free(a);
}

double my_roughenergy ( int * i1, double * x1, double * y1, double * z1, int n1, int * i2, 
                        double * x2, double * y2, double * z2, int n2, double cut,
                        double sigma, double epsilon ) {
   int i,j;
   double cut2=cut*cut;
   double d2,E=0.0;
   s6=sigma*sigma*sigma*sigma*sigma*sigma;
   rcut=pow(2,1./6.)*sigma
   if (rcut>cut) {
      rcut=cut
   }
   rcut2=rcut*rcut
   for (i=0;i<n1;i++) {
      for (j=0;j<n2;j++) {
         if (i1[i]!=i2[j]) {
           d2 =(x1[i]-x2[j])*(x1[i]-x2[j]);
           d2+=(y1[i]-y2[j])*(y1[i]-y2[j]);
           d2+=(z1[i]-z2[j])*(z1[i]-z2[j]);
           
           if (d2<rcut2) {
              di6=s6/(d2*d2*d2)
              di12=di6**2
              E+=4*epsilon*(di12-di6)+epsilon
           }
         }
      }
   }
   return E;
}
