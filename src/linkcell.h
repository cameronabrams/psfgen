#ifndef _LINKCELL_H_
#define _LINKCELL_H_

typedef struct LINKCELL {
    
    double xmin,xmax,xspan;
    double ymin,ymax,yspan;
    double zmin,zmax,zspan;
    double cut; // cutoff
    double dx, dy, dz; // cell dimensions in x, y , z
    int xnc, ync, znc; // number of cells along x, y, z
    int nc; // number of cells

    double * x, * y, * z; // particle position arrays, dim[n]
    int n; // number of particles
    int * cx, * cy, * cz; // cell coordinates of each particle, dim[n]
    int * ci; // linear cell index of each particle, dim[n]
    
    int **** pa; // particle-arrays dim[xnc][ync][znc][np[][][]]
    int *** np; // particle counts 
} linkcell;

linkcell * linkcell_new ( double * x, double * y, double * z, int n, double cut, int verbose );
linkcell * my_roughenergy_setup ( double * x2, double * y2, double * z2, int n2, double cellsize );
void linkcell_free ( linkcell * ls );
#endif
