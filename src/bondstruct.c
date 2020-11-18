#include "bondstruct.h"


bondstruct * new_bondstruct ( int * ia, int na, int nb ) {
   int i,j;
   bondstruct * bs = (bondstruct*)malloc(sizeof(bondstruct));
   //bs->nr=0;
   //bs->sr=0;
   bs->na=na;
   bs->mb=5; // by assumption
   // populate the list of atom-indexes
   bs->ia=(int*)malloc(na*sizeof(int));
   for (i=0;i<na;i++) bs->ia[i]=ia[i];
   //bs->rl=(int*)malloc(na*sizeof(int));
   //for (i=0;i<na;i++) bs->rl[i]=0;

   // set up blank lists of bond partners for each atom
   bs->ba=(int**)malloc(na*sizeof(int*));
   for (i=0;i<na;i++) {
     bs->ba[i]=(int*)malloc(bs->mb*sizeof(int));
     for (j=0;j<bs->mb;j++) bs->ba[i][j]=-1;
   }

   // set up blank lists of bonds
   bs->nb=nb;
   bs->isactive=(int*)malloc(nb*sizeof(int));
   bs->isrotatable=(int*)malloc(nb*sizeof(int));
   bs->bctr=0;
   bs->b=(int**)malloc(nb*sizeof(int*));
   for (i=0;i<nb;i++) {
      bs->b[i]=(int*)malloc(2*sizeof(int));
   }
   bs->bra=(int**)malloc(nb*sizeof(int*));
   for (i=0;i<nb;i++) {
      bs->bra[i]=(int*)malloc(na*sizeof(int));
      for (j=0;j<na;j++) bs->bra[i][j]=-1;
   }
   bs->bran=(int*)malloc(nb*sizeof(int));

   return bs;
}

void free_bondstruct ( bondstruct * bs ) {
   int i;
   if (bs->ia) free(bs->ia);
   //if (bs->rl) free(bs->rl);
   if (bs->ba) {
      for (i=0;i<bs->na;i++) free(bs->ba[i]);
      free(bs->ba);
   }
   if (bs->b) {
      for (i=0;i<bs->nb;i++) free(bs->b[i]);
      free(bs->b);
   }
   if (bs->bra) {
      for (i=0;i<bs->nb;i++) free(bs->bra[i]);
      free(bs->bra);
   }
   if (bs->bran) free(bs->bran);
   free(bs);
}

// get the index of atom with atom-index a in the bondstruct atom-index array
int bondstruct_getlocalatomindex ( bondstruct * bs, int a ) {
   int la;
   for (la=0;la<bs->na && bs->ia[la]!=a;la++);
   if (la==bs->na) {
     //printf("ERROR: atom-index %i is not in the bondstruct.\n",a);
     return -1;
   }
   return la;
}

void bondstruct_print ( bondstruct * bs ) {
   int i,j;
   printf("CFABOND/C) Array of rotatable bonds ('+' denotes active):\n");
   if (bs->b && bs->nb) {
      for (i=0;i<bs->nb;i++) {
         if (bs->isrotatable[i]) {
            printf("[%i] %c%i-%i : rs(%i)[",i,(bs->isactive[i]?'+':'-'),bs->b[i][0],bs->b[i][1],bs->bran[i]);
            for (j=0;j<bs->bran[i];j++) printf("%i ",bs->bra[i][j]);
            printf("]\n");
         }
      }
   }
/*    printf("CFABOND/C) Atom bondlists:\n");
   for (i=0;i<bs->na;i++) {
      printf("%i : ",bs->ia[i]);
      for (j=0;j<bs->mb&&bs->ba[i][j]!=-1;j++) printf("%i ",bs->ba[i][j]);
      printf("\n");
   }
 */
 }


int bondstruct_getbondindex ( bondstruct * bs, int i, int j ) {
   int k=0;
   int m=-1;
   for (k=0;k<bs->nb;k++) {
     if (bs->b[k][0]==i&&bs->b[k][1]==j) {
        m=k;
        break;
     }
   }
   return m;
}

int * bondstruct_getrightside_pointer ( bondstruct * bs, int b ) {
   return bs->bra[b];
}

int bondstruct_getrightside_count ( bondstruct * bs, int b ) {
   return bs->bran[b];
}

void bondstruct_addbond( bondstruct * bs, int i, int j) {
   if (bs->bctr==bs->nb) {
      printf("Error: too many bonds! (bug)\n");
      return;
   }
   //printf("assigning %i to 0th element of b[%i]\n",i,bs->bctr);fflush(stdout);
   bs->b[bs->bctr][0]=i;
   //printf("assigning %i to 1th element of b[%i]\n",j,bs->bctr);fflush(stdout);
   bs->b[bs->bctr][1]=j;
   //printf("setting rotatable flag\n");fflush(stdout);
   bs->isrotatable[bs->bctr]=1; // by assumption
   //printf("setting active flag\n");fflush(stdout);
   bs->isactive[bs->bctr]=1; // by assumption
   bs->bctr++;
}
// import bonds from TcL/vmd [$sel get index] and [$set getbonds] 
// a is an atom index and bl[] is an array of atom indices of its bonding partners
// as determined by [$sel getbonds].  nb is the number of such partners.
void bondstruct_importbonds ( bondstruct * bs, int a, int * bl, int nb ) {
  if (bs) {
    int i,ia;
    ia=bondstruct_getlocalatomindex(bs,a);
    if (nb>bs->mb) {
       printf("ERROR: atom %i has too many bonds\n",a);
       return;
    } else {
      for (i=0;i<nb;i++) {
         if (bl[i]==0) {
            printf("ERROR: position %i of atom %i neighbor list is zero!\n",i,a);
         }
         bs->ba[ia][i]=bl[i];
         //printf("adding bond %i %i...\n",a,bl[i]);fflush(stdout);
         bondstruct_addbond(bs,a,bl[i]);
      }
    }
  }
}

void bondstruct_setbond_rotatable ( bondstruct * bs, int a, int b, int flag ) {
   int k=bondstruct_getbondindex(bs,a,b);
   bs->isrotatable[k]=flag;
}

void bondstruct_setbond_active ( bondstruct * bs, int a, int b, int flag ) {
   int k=bondstruct_getbondindex(bs,a,b);
   bs->isactive[k]=flag;
}

// generate count of rotatable bonds and a map of index [0:nrb-1] to [0:nb-1]
void bondstruct_maprotatables ( bondstruct * bs ) {
   int i,j;
   bs->nrb=0;
   for (i=0;i<bs->nb;i++) bs->nrb+=bs->isrotatable[i]?1:0;
   bs->r2b=(int*)malloc(bs->nrb*sizeof(int));
   j=0;
   for (i=0;i<bs->nb;i++) {
      if (bs->isrotatable[i]) {
         bs->r2b[j]=i;
         j++;
      }
   }
}

int bondstruct_r2b ( bondstruct * bs, int r ) {
   return bs->r2b[r];
}

int bondstruct_getnb ( bondstruct * bs ) {
   return bs->nb;
}

int bondstruct_getrnb ( bondstruct * bs ) {
   return bs->nrb;
}

int * bondstruct_getbondpointer ( bondstruct * bs, int i ) {
   return bs->b[i];
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
   
/*    // allocate the array of arrays
   bs->bra=(int**)malloc(bs->nb*sizeof(int*));
   for (i=0;i<bs->nb;i++) {
     bs->bra[i]=(int*)malloc(bs->na*sizeof(int));
     for (j=0;j<bs->na;j++) bs->bra[i][j]=-1;
   }
   // allocate array of right-side counts
   bs->bran=(int*)malloc(bs->nb*sizeof(int));
 */   // for each bond, identify and save all right-side atoms
//   printf("#### in makerightsides processing %i bonds\n",bs->nb);fflush(stdout);

   // for each bond
   for (k=0;k<bs->nb;k++) {
      // get indices of two atoms; 'a' is the left-atom index and 'b' is the right-atom index
      a=bs->b[k][0];
      b=bs->b[k][1];
   //   printf("#### in makerightsides at bond %i : %i %i\n",k,a,b);fflush(stdout);
      // get the two local indices
      la=bondstruct_getlocalatomindex(bs,a);
      lb=bondstruct_getlocalatomindex(bs,b);
   //   printf("#### local indices %i %i\n",la,lb);fflush(stdout);
      // put all atoms b is bonded to on the right-side list _except_ atom a
      bs->bran[k]=0;
      if (!bs->isrotatable[k]) continue;
      for (i=0;i<bs->mb && bs->ba[lb][i]!=-1;i++) {
         if (bs->ba[lb][i] != a) {
   //         printf("#### added %i-neighbor %i to bond-%i rightside list\n",b,bs->ba[lb][i],k);fflush(stdout);
            bs->bra[k][bs->bran[k]++]=bs->ba[lb][i];
         }
      }
     // enter the grow-loop; for each element on bs->bra[k][], add _its_ bond partners to this
     // right-side array, provided they are not already in it
   //   printf("#### growing rightside of bond %i\n",k);fflush(stdout);
      grow=1;
      while (grow) {
         grow=0;
         // save the current right-side array size
         // lnr=bs->bran[k];
         // for each atom on this right-side array
         for (i=0;i<bs->bran[k];i++) {
            // get the local index
            li=bondstruct_getlocalatomindex(bs,bs->bra[k][i]);
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

int bondstruct_isactive ( bondstruct * bs, int b ) {
   return bs->isactive[b] && bs->bran[b]>0;
}

// deactivate bonds whose right-side lists contain the index 'fa'
void bondstruct_deactivate_by_fixed ( bondstruct * bs, int fa ) {
   int i,j;
   for (i=0;i<bs->nb;i++) {
      for (j=0;j<bs->bran[i]&&bs->bra[i][j]!=fa;j++);
      if (j<bs->bran[i]) bs->isactive[i]=0;
   }
}

int bondstruct_arebonded ( bondstruct * bs, int a, int b ) {
   int * ba;
   int i;
   int la=bondstruct_getlocalatomindex(bs,a);
   if (la==-1) return 0;
   ba=bs->ba[la];
   for (i=0;i<bs->mb&&ba[i]!=-1&&ba[i]!=b;i++);
   if (i==bs->mb) return 0;
   if (ba[i]==-1) return 0;
   return 1;
}
 
void free_intarray ( int * a ) {
   free(a);
}

int isin ( int * arr, int n, int t ) {
   int i;
   for (i=0;i<n&&arr[i]!=t;i++);
   if (i<n) return 1;
   return 0;
}


double my_roughenergy ( int * i1, double * x1, double * y1, double * z1, int n1, 
                        int * i2, int n2,
                        double cut, double sigma, double epsilon, bondstruct * bs, linkcell * ls ) {
   int i,j,ij;
   int a,b,icx,icy,icz,tx,ty,tz,dx,dy,dz,*pa,np;
   double d2,E=0.0,di6,di12,bx,by,bz;
   double s6=sigma*sigma*sigma*sigma*sigma*sigma;
   double rcut=pow(2,1./6.)*sigma,rcut2;
   if (rcut>cut) {
      rcut=cut;
   }
   rcut2=rcut*rcut;
   for (i=0;i<n1;i++) {
      // pairwise sum over all sel1 pairs
      a=i1[i];
      for (ij=i+1;ij<n1;ij++) {
         if (bondstruct_arebonded(bs,a,a)) continue;
         d2 =(x1[i]-x1[ij])*(x1[i]-x1[ij]);
         d2+=(y1[i]-y1[ij])*(y1[i]-y1[ij]);
         d2+=(z1[i]-z1[ij])*(z1[i]-z1[ij]);
         if (d2<rcut2) {
            di6=s6/(d2*d2*d2);
            di12=di6*di6;
            E+=4*(di12-di6)+1;
         }      
      }
      // linkcell sum over all sel1-sel2 pairs
      icx=(int)((x1[i]-ls->xmin)/ls->dx);
      icy=(int)((y1[i]-ls->ymin)/ls->dy);
      icz=(int)((z1[i]-ls->zmin)/ls->dz);
      for (dx=-1;dx<2;dx++) {
         tx=icx+dx;
         if (tx==ls->xnc||tx==-1) continue;
         for (dy=-1;dy<2;dy++) {
            ty=icy+dx;
            if (ty==ls->ync||ty==-1) continue;
            for (dz=-1;dz<2;dz++) {
               tz=icy+dz;
               if (tz==ls->znc||tz==-1) continue;
               pa=ls->pa[tx][ty][tz];
               np=ls->np[tx][ty][tz];
               for (j=0;j<np;j++) {
                  b=i2[pa[j]];
                  if (bondstruct_arebonded(bs,a,b)) continue;
                  bx=ls->x[pa[j]];
                  by=ls->y[pa[j]];
                  bz=ls->z[pa[j]];
                  d2 =(x1[i]-bx)*(x1[i]-bx);
                  d2+=(y1[i]-by)*(y1[i]-by);
                  d2+=(z1[i]-bz)*(z1[i]-bz);
                  if (d2<rcut2) {
                     di6=s6/(d2*d2*d2);
                     di12=di6*di6;
                     E+=4*(di12-di6)+1;
                  }
               }
            }
         }
      }
   }
   return E*epsilon;
}

double max ( double * a, int n ) {
   double r = -1.e23;
   int i;
   for (i=0;i<n;i++) {
      if (a[i]>r) {
         r=a[i];
      }
   }
   return r;
}

double min ( double * a, int n ) {
   double r = 1.e23;
   int i;
   for (i=0;i<n;i++) {
      if (a[i]<r) {
         r=a[i];
      }
   }
   return r;
}

linkcell * linkcell_new ( double * x, double * y, double * z, int n, double cut ) {
    int i, j, k;
    int icx, icy, icz, pi, lci;
    int * npl;

    linkcell * ls=(linkcell*)malloc(sizeof(linkcell));
    ls->n;
    ls->x=x;
    ls->y=y;
    ls->z=z;
    ls->ci=(int*)malloc(n*sizeof(int));
    ls->cx=(int*)malloc(n*sizeof(int));
    ls->cy=(int*)malloc(n*sizeof(int));
    ls->cz=(int*)malloc(n*sizeof(int));
    
    ls->xmin=min(x,n)-cut;
    ls->xmax=max(x,n)+cut;
    ls->xspan=ls->xmax-ls->xmin;
    ls->xnc=(int)(ls->xspan/ls->cut+1);
    ls->dx=ls->xspan/ls->xnc;

    ls->ymin=min(y,n)-cut;
    ls->ymax=max(y,n)+cut;
    ls->yspan=ls->ymax-ls->ymin;
    ls->ync=(int)(ls->yspan/ls->cut+1);
    ls->dy=ls->yspan/ls->ync;

    ls->zmin=min(z,n)-cut;
    ls->zmax=max(z,n)+cut;
    ls->zspan=ls->zmax-ls->zmin;
    ls->znc=(int)(ls->zspan/ls->cut+1);
    ls->dz=ls->zspan/ls->znc;

    ls->nc=ls->xnc*ls->ync*ls->znc;

    printf("Sel2 linkcell %i particles with cut %.4f:\n",ls->n,ls->cut);
    printf("   x: [%.4f,%.4f,%.4f] %i\n",ls->xmin,ls->dx,ls->xmax,ls->xnc);
    printf("   y: [%.4f,%.4f,%.4f] %i\n",ls->ymin,ls->dy,ls->ymax,ls->ync);
    printf("   z: [%.4f,%.4f,%.4f] %i\n",ls->zmin,ls->dz,ls->zmax,ls->znc);
    printf("   total cells: %i\n",ls->nc);

    npl=(int*)malloc(ls->nc*sizeof(int));
    for (i=0;i<ls->nc;i++) {
        npl[i]=0;
    }

    for (i=0;i<n;i++) {
        ls->cx[i]=(int)((ls->x[i]-ls->xmin)/ls->dx);
        ls->cy[i]=(int)((ls->y[i]-ls->ymin)/ls->dy);
        ls->cz[i]=(int)((ls->z[i]-ls->zmin)/ls->dz);
        ls->ci[i]=ls->cx[i]*ls->ync*ls->znc+ls->cy[i]*ls->znc+ls->cz[i];
        npl[ls->ci[i]]++;
    }

    ls->pa=(int****)malloc(ls->xnc*sizeof(int***));
    ls->np=(int***)malloc(ls->xnc*sizeof(int**));
    for (i=0;i<ls->xnc;i++) {
        ls->pa[i]=(int***)malloc(ls->ync*sizeof(int**));
        ls->np[i]=(int**)malloc(ls->ync*sizeof(int*));
        for (j=0;j<ls->ync;j++) {
            ls->pa[i][j]=(int**)malloc(ls->znc*sizeof(int*));
            ls->np[i][j]=(int*)malloc(ls->znc*sizeof(int));
            for (k=0;k<ls->znc;k++) {
                lci=i*ls->ync*ls->znc+j*ls->znc+k;
                ls->pa[i][j][k]=(int*)malloc(npl[lci]*sizeof(int*));
                ls->np[i][j][j]=0;
            }
        }
    }
    for (i=0;i<n;i++) {
        icx=ls->cx[i];
        icy=ls->cy[i];
        icz=ls->cz[i];
        pi=ls->np[icx][icy][icz];
        ls->pa[icx][icy][icz][pi]=i;
        ls->np[icx][icy][icz]++;
    }
    return ls;
}

linkcell * my_roughenergy_setup ( double * x2, double * y2, double * z2, int n2, double cut ) {
   linkcell * ls = linkcell_new(x2,y2,z2,n2,cut);
   return ls;
}