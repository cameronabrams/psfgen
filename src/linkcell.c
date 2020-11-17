#include "linkcell.h"

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

    linkcell * ls=malloc(sizeof(linkcell));
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