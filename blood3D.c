/* Bo-Göran Wallner */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "nr3.h"
#include "transforms.c"
#include "memoryallocation2.c"


#define twopi 6.283185307177959
#define pi 3.14159265358

#define DEBUG 1
#define MAX_LIMIT 20000
#define PRINTOUT 0
#define FIELD_MAX_REAL 100000
#define FIELD_MAX_IM 100000

/* Choose angle setup for different simulations */
#define ANGLES_RANDOM         1
#define ANGLE_THETA_ZERO      0
#define ANGLE_THETA_PI_HALF   0
#define ANGLE_THETA_PI_FOURTH 0

typedef struct
{
	int amp, xc, yc, zc, dx, dy;
	double theta, fi, psi;
} bloodData;

typedef struct 
{
	int xAnt, yAnt, zAnt, iAnt, xout, zout;
	double dx, dy, dz, fraction;
	int xbox, ybox, zbox;
} sampData;


/* Describes medium */
typedef struct
{
	double afreq; /* angular frequency */
	double epsilonRe, epsilonIm; /* inclusion permittivity */
	double backRe, backIm; /* backgound permittivity */
} modelData;


/* output of structures ************************************/

void printfsampData(sampData *sD)
{
printf("Discretization data\n");
printf("xAnt = %d\t, yAnt = %d\t, zAnt = %d\n",sD->xAnt, sD->yAnt, sD->zAnt);
printf("dx = %.4f\t, dy = %.4f\t, dz = %.4f\n\n",sD->dx, sD->dy, sD->dz);
}


void fprintfsampData(FILE *fp, sampData *sD)
{
	fprintf(fp, "Discretization data\n");
	fprintf(fp, "xAnt = %d\t, yAnt = %d\t, zAnt = %d\n",
	sD->xAnt, sD->yAnt, sD->zAnt);
	fprintf(fp, "dx = %.4f\t, dy = %.4f\t, dz = %.4f\n\n",
	sD->dx, sD->dy, sD->dz);
}


void printfmodelData(modelData *mD)
{
	printf("Modelling data\n");
	printf("angular frequency = %.4f\n", mD->afreq);
	printf("inclusion = %.4f-i %.4f\n", mD->epsilonRe, mD->epsilonIm);
	printf("backgroud = %.4f-i %.4f\n", mD->backRe, mD->backIm);
}


void fprintfmodelData(FILE *fp, modelData *mD)
{
	fprintf(fp, "Modelling data\n");
	fprintf(fp, "angular frequency = %.4f\n", mD->afreq);
	fprintf(fp, "inclusion = %.4f-i %.4f\n", mD->epsilonRe, mD->epsilonIm);
	fprintf(fp, "backgroud = %.4f-i %.4f\n", mD->backRe, mD->backIm);
}

/* some complex valued algebra *****************************/
void compsqrt(double a, double b, double *c, double *d)
{
	*d = sqrt(a*a + b*b);
	*c = sqrt((*d + a)/2.0);
	*d = sqrt((*d -a )/2.0);
	*d = (b < 0) ? -*d : *d;
}

void compmult(double ar, double ai, double br, double bi, double *cr, double *ci)
{
	*cr = ar*br - ai*bi;
	*ci = ar*bi + ai*br;

    /* Assert if out of boundaries */
	assert(sizeof(*cr) == sizeof(double));
	assert(sizeof(*ci) == sizeof(double));
}

void compdiv(double ar, double ai, double br, double bi, double *cr, double *ci)
{
	*ci = br*br + bi*bi;
	*cr = (ar*br + ai*bi)/(*ci);
	*ci = (br*ai - bi*ar)/(*ci);
}

double compabs(double ar, double ai)
{
	return sqrt(ar*ar + ai*ai);
}

double absolut(double a)
{
	return (a < 0) ? -a : a;
}

double random1(double a){
	return a*rand()/(1.0*RAND_MAX);
}

/* 
 * Memoryallocation of bloodData 3D array
 */
bloodData ****mallocBloodData3DArray(int dim1, int dim2, int dim3)
{
	int i, j, k;
	bloodData ****a = (bloodData ****) malloc(dim1 * sizeof(bloodData ***));
	for (i=0; i<dim1; i++){
		a[i] = (bloodData ***) malloc(dim2 * sizeof(bloodData **));
		for (j=0; j<dim2; j++){
			a[i][j] = (bloodData **) malloc(dim3 * sizeof(bloodData *));
			for (k=0; k<dim3; k++) {
					a[i][j][k] = (bloodData *) malloc(sizeof(bloodData));
			}
		}
	}
	return a;
}

/* 
 * Create 3D array with bloodData
 */
void initBloodData3DArray(sampData *sD, bloodData ****bD)
{
	int nbrOfXbox, nbrOfYbox, nbrOfZbox;
	int zb, xb, yb;

    /* Calculate number of cells. E.g if we use 1024 sample points    */
	/* and each cell is 128 points we get 1024/128 = 8 cells          */
	nbrOfXbox = sD->xAnt/sD->xbox;
	nbrOfZbox = sD->zAnt/sD->zbox;
	nbrOfYbox = sD->yAnt/sD->ybox;

	for (zb = 0; zb < nbrOfZbox; zb++)
	{
		/* Displace centers of whole xy-plane layer */
		/* to create randomness in geometry.        */
		bD[zb][yb][xb]->dx = (int)random1(sD->xbox);
		bD[zb][yb][xb]->dy = (int)random1(sD->ybox);

		for (yb = 0; yb < nbrOfYbox; yb++)
		{
			for (xb = 0; xb < nbrOfXbox; xb++)
			{
				bD[zb][yb][xb]->zc = (zb+0.5)*sD->zbox;
				bD[zb][yb][xb]->xc = (xb+0.5)*sD->xbox;
				bD[zb][yb][xb]->yc = (yb+0.5)*sD->ybox;
				bD[zb][yb][xb]->amp = 1;

				/* E.g. if we have 1024 samplepoints and choose 128 size cells        */
				/* we get 8 cells having centers (0+0.5)*128 = 64, (1+0.5)*128 = 192, */
				/* (2+0.5)* 128 = 320, ... (7+0.5)*128 = 960                          */

                /* Assign values for the Euler angles */
#ifdef (ANGLES_RANDOM == 1)
				bD[zb][yb][xb]->theta = random1(2*pi);
				bD[zb][yb][xb]->fi = random1(2*pi);
				bD[zb][yb][xb]->psi = random1(2*pi);
#endif

#ifdef (ANGLE_THETA_ZERO == 1)
                bD[zb][yb][xb]->theta = 0;
				bD[zb][yb][xb]->fi = random1(2*pi);
				bD[zb][yb][xb]->psi = random1(2*pi);
#endif

#ifdef (ANGLE_THETA_PI_HALF == 1)
                bD[zb][yb][xb]->theta = pi/2;
				bD[zb][yb][xb]->fi = random1(2*pi);
				bD[zb][yb][xb]->psi = random1(2*pi);
#endif

#ifdef (ANGLE_THETA_PI_FOURTH == 1)
                bD[zb][yb][xb]->theta = pi/4;
				bD[zb][yb][xb]->fi = random1(2*pi);
				bD[zb][yb][xb]->psi = random1(2*pi);
#endif
			}
		}
	}
}

/* 
* Fill 2D geometry layer fox fixed z 
*/

static void create2DGeometry(sampData *sD, bloodData ***bD, MatInt& geometry2D, int z)
{
	/* Coordinates in global coordinate system */
	int x, y, zb, xb, dx, dy, yb;
	/* Coordinates in local cell with disk in origo, e.g. xl = f(x,y,z) */
	int xl, yl, zl;
	/* Coordinates in Euler rotated system, e.g. xe = f(xl,yl,zl)       */
	int xe, ye, ze, Re;
    /* Blood-cell constants */
	double C0, C2, C4, R0, D_xDiff;
	/* Trigonometric values */
	double cost, sint, sinf, cosf, sinp, cosp;

    /* Iterate over the whole xy-plane */
	for (y=0; y<sD->yAnt; y++) 
	{
		for (x=0; x<sD->xAnt; x++) 
		{
			/* Calculate which cell */
			zb=z/sD->zbox;
			xb=x/sD->xbox;
			yb=y/sD->ybox;

			/* Pre-assign Euler angels to simplify expressions */
			cost = cos(bD[yb][xb]->theta);
			sint = sin(bD[yb][xb]->theta);
			sinf = sin(bD[yb][xb]->fi);
			cosf = cos(bD[yb][xb]->fi);
			cosp = cos(bD[yb][xb]->psi);
			sinp = sin(bD[yb][xb]->psi);

			/* The task is now to determine if a point (x,y) is within  */
            /* or outside of the cross-section of the disk for a given  */
			/* fixed z. To do this we need to transform the (x,y,z)     */
			/* global "entire geometrical model" coordinates into a     */
			/* local axial-symmetric coordinate systems for which we    */
			/* know the equation of the disk.                           */

            /* Start with transformation from global to local coordinate */
			/* system of the center of the disk in the cell (xc, yc) by  */
			/* (xl,yl,zl) = (x,y,z) - (xc,yc,zc)                         */
			xl = x - bD[yb][xb]->xc;
			yl = y - bD[yb][xb]->yc;
			zl = z - bD[yb][xb]->zc;

			/* We are now at disk origo. Apply Euler rotation to rotate */
			/* coordinate system into the correct angles and describe   */
			/* the disk in local coordinates (xe,ye,ze).                */ 
			xe = xl*(cosp*cosf-cost*sinf*sinp) + yl*(cosp*sinf+cost*cosf*sinp) + zl*(sinp*sint);
			ye = xl*(-sinp*cosf-cost*sinf*cosp) + yl*(-sinp*sinf+cost*cosf*cosp) + zl*(cosp*sint);
			ze = xl*(sint*sinf) + yl *(-sint*cosf) + zl*(cost);

            /* So far we have not assumed anything about the exact     */
			/* geometry. We have barely stated, we have any object     */
			/* shape in (x,y,z) and translate this object into a       */
			/* local coordinate system (xe,ye,ze) by translation       */
			/* and Euler rotation (since we decided to describe the    */
			/* object with Euler angles).                              */

            /************ BICONCAVE DISK / BLOODCELL APPROACH **********/

            /* Suggestion by A.Karlsson "Numerical Simulations of Light Scattering          */
			/* by Red Blood Cells" the equation for the biconcave disk can be               */
			/* described as: D(x) = Sqrt(1-(x/R0)^2) * (C0+C2*(x/R0)^2+C4*(x/R0)^4))        */
			/* where C0=0.81um, C2=7.83um, C4=-4.39um and R0=3.91um and this corresponds    */
			/* to a biconcave disk with diameter a=7.76um and max thickness b=2.55um.       */
			/* Thus, the a=7.76um shall correspond to sD->xbox samplepoints.                */
			/* We apply the following tranformations: x -> Sqrt(xe^2+ye^2) and for D(x)     */
			/* we set D(x)^2 -> (2ze)^2 = 4*ze^2. This give a 3D expression and we can use  */
			/* 4*ze^2 - 4*(1-(Re/R0)^2) * (C0+C2*(Re/R0)^2+C4*(Re/R0)^4))^2 > 0 and if      */
			/* this expression holds then the ze is *not* within the biconcave disk.        */
			/* Also Re=sqrt(xe^2+ye^2) < half cell size to determine if (xe,ye) is within   */
			/* the disk or not.                                                             */
			R0 = 64.4948;
			C0 = 13.3608;
			C2 = 129.1546;
			C4 = -72.41243;
            /* Calculate (xe,ye) corresponding radius */
			Re = sqrt(xe^2 + ye^2);

			/* The whole xy-layer shall be translated */
			dx = bD[yb][xb]->dx;
			dy = bD[yb][xb]->dy;

			/* Check if (xe,ye) is within disk since expression only valid within */
			/* it and that exclude e.g. corners of the cell.                      */
			if(Re < sD->xbox/2)
			{
				/* Re is valid to use in the expression */
				D_xDiff = ze^2 - (1-(Re/R0)^2) * (C0 + C2*(Re/R0)^2 + C4 * (Re/R0)^4)^2;
				/* Check if ze is outside or within the biconcave disk */
				if ( D_xDiff > 0.0)
				{
					/* ze is outside the disk */
					geometry2D[(x+dx)%sD->xAnt][(y+dy)%sD->yAnt] = 0;
				}
				else
				{
					/* ze is within the biconcave disk */
					geometry2D[(x+dx)%sD->xAnt][(y+dy)%sD->yAnt] = 1;
				}
			}
			else
			{
				/* Re is outside the disk */
				geometry2D[(x+dx)%sD->xAnt][(y+dy)%sD->yAnt] = 0;
			}
		}
	}
}

/* migrate the fields one layer ****************************************/
void migrate(sampData *sD, modelData *mD, MatInt& slow, MatDoub& u1)
{
	int xmax=sD->xAnt, zmax=sD->zAnt, imax=sD->iAnt, ymax=sD->yAnt;
	double dx=sD->dx, dz=sD->dz, dy=sD->dy ;
	double w=mD->afreq, eta=0.0;
	double s, xiX, xiY, smin=1000.0, smax=0.0, b, temp, propr, propi;
	double ii;
	double epsr, epsi, ga2r, ga2i, gar, gai, kr, ki, k2r;
	int x, i, y, pr=0, ev=0;

	vector<MatDoub> v1(2);
	for (Int i=0; i<2;i++) v1[i].resize(ymax, 2*xmax);

    /***********/
//    for (y=0; y<ymax; y++) {
//        for (x=0; x<2*xmax; x+=2) {
//            printf("Re(umig(%d %d)) : %f \n", y, x, u1[y][x]);
//            printf("Im(umig(%d %d)) : %f \n", y, x+1, u1[y][x+1]);
//        }
//    }
    /*************/

	ddfourn(u1, xmax, ymax);
    
	for (y=0; y<ymax; y++) 
	{
		for (x=0; x<xmax; x++){ /* slowness max and min */
			temp = slow[y][x];
			smax = (smax > temp) ? smax : temp;
			smin = (smin < temp) ? smin : temp;
		}
	}

	for (i=0; i<imax; i++) 
	{   /*interpolation slowness */
		s = smin + i*(smax-smin)/(imax-1); /* slowness */
		epsr = 1.0454 + s*mD->epsilonRe; /* Re permittivity */
		epsi = s*mD->epsilonIm; /* Im permittivity*/
		compmult(eta, w, eta, w, &kr, &ki); /* wave number, c^{-1} s */
		compmult(kr, ki, epsr, epsi, &k2r, &ga2i); /* sqr wave number, .^{2} */
		xiX = 0.0; /* zero transverse wave number */
		xiY = 0.0;
        
        /* Must handle y=0 separately since (ymax-y) exist for y>0 and y=0..ymax-1 */
		y = 0.0;
		ga2r = k2r + xiX*xiX + xiY*xiY;
		compsqrt(ga2r, ga2i, &gar, &gai);
		propr = cos(-dz*gai) * exp(-dz*gar);
		propi = sin(-dz*gai) * exp(-dz*gar);
		compmult(propr, propi, u1[y][0], u1[y][1], &v1[i][y][0], &v1[i][y][1]);
        
        for (x=2; x<xmax+2; x+=2){
            /* zero transverse wave number */
            xiX = x/(xmax*2.0*dx)*twopi;
            ga2r = k2r + xiX*xiX + xiY*xiY;
            compsqrt(ga2r, ga2i, &gar, &gai);
            propr = cos(-dz*gai) * exp(-dz*gar);
            propi = sin(-dz*gai) * exp(-dz*gar);
            compmult(propr, propi, u1[y][x], u1[y][x+1], &v1[i][y][x], &v1[i][y][x+1]);
            compmult(propr, propi, u1[y][2*xmax-x], u1[y][2*xmax-x+1],&v1[i][y][2*xmax-x],
                     &v1[i][y][2*xmax-x+1]);
        }

        /* And now handle y>0 */
		for (y=1; y<ymax/2-1; y++)
		{
			xiX = 0.0; /* zero transverse wave number */
			xiY = y/(ymax*dy)*twopi;
			ga2r = k2r + xiX*xiX + xiY*xiY;
			compsqrt(ga2r, ga2i, &gar, &gai);
			propr = cos(-dz*gai) * exp(-dz*gar);
			propi = sin(-dz*gai) * exp(-dz*gar);
			compmult(propr, propi, u1[y][0], u1[y][1], &v1[i][y][0], &v1[i][y][1]);
		    compmult(propr, propi, u1[ymax-y][0], u1[ymax-y][1], &v1[i][ymax-y][0],
			&v1[i][ymax-y][1]);
			for (x=2; x<xmax+2; x+=2)
			{
				xiX = x/(xmax*2.0*dx)*twopi;
				ga2r = k2r + xiX*xiX + xiY*xiY;
				compsqrt(ga2r, ga2i, &gar, &gai);
				propr = cos(-dz*gai) * exp(-dz*gar);
				propi = sin(-dz*gai) * exp(-dz*gar);
				compmult(propr, propi, u1[y][x], u1[y][x+1], &v1[i][y][x], &v1[i][y][x+1]);
				compmult(propr, propi, u1[y][2*xmax-x], u1[y][2*xmax-x+1],
				&v1[i][y][2*xmax-x], &v1[i][y][2*xmax-x+1]);
				compmult(propr, propi, u1[ymax-y][x], u1[ymax-y][x+1],
				&v1[i][ymax-y][x], &v1[i][ymax-y][x+1]);
				compmult(propr, propi, u1[ymax-y][2*xmax-x], u1[ymax-y][2*xmax-x+1],
				&v1[i][ymax-y][2*xmax-x],
				&v1[i][ymax-y][2*xmax-x+1]);
			}
		}

		/* Make inverse FFT */
		ddfourninv(v1[i], xmax, ymax);
					    /***********/
//    for (y=0; y<ymax; y++) {
//        for (x=0; x<2*xmax; x+=2) {
//            printf("Re(umig(%d %d)) : %f \n", y, x, v1[i][y][x]);
//            printf("Im(umig(%d %d)) : %f \n", y, x+1, v1[i][y][x+1]);
//        }
//    }
    /*************/
        
    } /* END for (i=0; i<imax; i++) { /*interpolation slowness */

    if ((smax-smin>0.001) && (imax>1))
    {
        for (y=0; y<ymax; y++) {
            for (x=0; x<xmax; x++) { /* interpolation */
                b = modf((slow[y][x]-smin)/(smax-smin)*(imax-1), &ii);
                i = ii;
                if (i<imax-1){
                    u1[y][2*x] = ((1.0-b)*v1[i][y][2*x] + b*v1[i+1][y][2*x]);
                    u1[y][2*x+1] = ((1.0-b)*v1[i][y][2*x+1] + b*v1[i+1][y][2*x+1]);
                }
                else {
                    u1[y][2*x] = v1[i][y][2*x];
                    u1[y][2*x+1] = v1[i][y][2*x+1];
                }
            }
        }
    }
    else {
        for (y=0; y<ymax; y++) {
            for (x=0; x<2*xmax; x++){ /*interpolation */
                u1[y][x] = v1[0][y][x];
            }
        }
    }

	/* Deallocate temp */
	//multidealloc3(v1, imax, ymax, 2*xmax);	
}


void interaction(sampData *sD, modelData *mD, MatInt &slow1, MatInt &slow2,
MatDoub& uin, MatDoub& urefl) {

	int xmax=sD->xAnt, zmax=sD->zAnt, ymax=sD->yAnt;
	double dx=sD->dx, dz=sD->dz, dy=sD->dy;
	double w=mD->afreq, eta=0.0;
	double epsr, epsi, kr, ki, k2r, gar, gai, ga2r, ga2i;
	double bsr, bsi, bkr, bki, bk2r, bgar, bgai, bga2r, bga2i;
	double refr, refi, xiX, xiY;
	int x, y, **refp, refl=0;
	
	refp = multialloc2int(ymax,xmax);
	for (y=0; y<ymax; y++) {
		for (x=0; x<xmax; x++){
			refp[y][x] = 0;
			if (slow1[y][x]-slow2[y][x] > 0.5){
				refp[y][x] = 1;
				refl = 1;
			}
			else if (slow2[y][x]-slow1[y][x] > 0.5){
				refp[y][x] = -1;
				refl = 1;
			}
		}
	}

	/* Obs: l�ngd e nu "ymax x 2*xmax" */
	if (refl == 1) {
		
		y = 0;
		
		/* Transform in-field */
		ddfourn(uin, xmax, ymax); /* spatial FFT */
		
		/* Set permitivity */
		epsr = 1.0454 + mD->epsilonRe; /* Re permittivity */
		epsi = mD->epsilonIm; /* Im permittivity*/

		/* Compute wavenumber for cells */
		compmult(eta, w, eta, w, &kr, &ki); /* wavenumber, c^{-1} s */
		compmult(kr, ki, epsr, epsi, &k2r, &ga2i); /* sqr wavenumber, .^{2} */
		
		/* Include background */
		bsr = 1.0; /* Re of the background */
		bsi = 0.0;

		/* Compute wavenumber for background */
		compmult(bsr, bsi, eta, w, &bkr, &bki); /* wavenumber, c^{-1} s */
		compmult(bkr, bki, bkr, bki, &bk2r, &bga2i); /* sqr wavenumber, .^{2} */
		
		xiY = 0.0;
		xiX = 0.0; /* zero transverse wave number */
		ga2r = k2r + xiX*xiX + xiY*xiY;
		compsqrt(ga2r, ga2i, &gar, &gai);
		bga2r = bk2r + xiX*xiX + xiY*xiY;
		compsqrt(bga2r, bga2i, &bgar, &bgai);
		compdiv(gar-bgar, gai-bgai, gar+bgar, gai+bgai, &refr, &refi);
		compmult(refr, refi, uin[y][0], uin[y][1], &urefl[y][0],
		&urefl[y][1]);

		for (y=0; y<ymax; y++) {
			for (x=0; x<xmax+2; x+=2){
				xiX = x/(xmax*2.0*dx)*twopi;
				xiY = y/(ymax*dy)*twopi;
				ga2r = k2r + xiX*xiX + xiY*xiY;
				compsqrt(ga2r, ga2i, &gar, &gai);
				bga2r = bk2r + xiX*xiX + xiY*xiY;
				compsqrt(bga2r, bga2i, &bgar, &bgai);
				compdiv(gar-bgar, gai-bgai, gar+bgar, gai+bgai, &refr, &refi);
				compmult(refr, refi, uin[y][x], uin[y][x+1], &urefl[y][x],
				&urefl[y][x+1]); /* reflected field, ur = R*uin */
				compmult(refr, refi, uin[y][2*xmax-x], uin[y][2*xmax-x+1],
				&urefl[y][2*xmax-x], &urefl[y][2*xmax-x+1]);
			}
		}

		ddfourninv(urefl, xmax, ymax); /* inverse spatial FFT */
		ddfourninv(uin, xmax, ymax); /* inverse spatial FFT fix ????? */
		for (y=0; y<ymax; y++){
			for (x=0; x<xmax; x++){
				urefl[y][2*x] = refp[y][x]*urefl[y][2*x];
				urefl[y][2*x+1] = refp[y][x]*urefl[y][2*x+1];
				uin[y][2*x] = uin[y][2*x] + urefl[y][2*x];
				uin[y][2*x+1] = uin[y][2*x+1] + urefl[y][2*x+1];
			}
		}
	}
	else {
		for (y=0; y<ymax; y++){
				for (x=0; x<xmax; x++){
				urefl[y][2*x] = 0.0; /* transmitted field, ut=T*uin = (1+R) uin */
				urefl[y][2*x+1] = 0.0;
			}
		}
	}

	/* Deallocate */
	multidealloc2int(refp, ymax, xmax);
}


/* main program to propagate the field through the a sequence of layers */
void propagate(sampData *sD, modelData *mD) {

	int x, y, z;
	double powerTransmitted;

	/* Indexing [z][y][x] will be used where z is in propagation direction  */
	/* Fields in plane [y][x] will have Re[field] = field[y][x]             */
	/* and Im[field] = field[y][x+1]                                        */

   /* Declare geometry layers */
    MatInt sampleLayer1(sD->yAnt, sD->xAnt);
    MatInt sampleLayer2(sD->yAnt, sD->xAnt);
    
	/* Declare field, Real + Imaginary parts */
    MatDoub umig(sD->yAnt, 2*sD->xAnt);
    MatDoub uref(sD->yAnt, 2*sD->xAnt);

    /* bD pointer to 3D array with bloodData */
	bloodData ****bD;

	/* Allocate memory for 3D array, every point containing bloodData */
	bD = mallocBloodData3DArray( sD->zAnt/sD->zbox, sD->yAnt/sD->ybox, sD->xAnt/sD->xbox );

	/* Initialize the 3D grid, blood cell centers only */
	initBloodData3DArray( sD, bD );

    /* Assign field init values.            */
	/* Input field:     Re = 1.0 Im = 0.0   */
	/* Reflected field: Re = 0.0 Im = 0.0   */
	for (y=0; y<sD->yAnt; y++)
	{
		for (x=0; x<2*sD->xAnt; x+=2)
		{
			umig[y][x] =   1.0;
			umig[y][x+1] = 0.0;
			uref[y][x]   = 0.0;
			uref[y][x+1] = 0.0;
		}
	}

	/* Start propagate from z=0 to z=(sD->zAnt-2) */
	for (z=0; z<(sD->zAnt - 1); z++)
	{
		/* Fill two layers in xy-plane with sample points */ 
		create2DGeometry(sD, bD[z/sD->zbox], sampleLayer1, z);
		create2DGeometry(sD, bD[(z+1)/sD->zbox], sampleLayer2, z+1);

		/* Migrate between the two layers */
		migrate(sD, mD, sampleLayer1, umig);
        
		/* Calculate interaction between two layers */
		interaction(sD, mD, sampleLayer1, sampleLayer2, umig, uref);

		/* Calculate power in z by summing over x-y plane */
		/* (uRe^2+uIm^2) * (1+epsilon[y][x/2]) */  /* TODO: Why x/2 ??? */
		powerTransmitted = 0.0;
		for ( y = 0; y < sD->yAnt; y++) 
		{
			for (x = 0; x < 2*sD->xAnt; x+=2) 
			{
				powerTransmitted += (umig[y][x]*umig[y][x]+umig[y][x+1]*umig[y][x+1])*
				                    (1.0+mD->epsilonRe*sampleLayer1[y][x/2]);
			}
		}

		/* normalize */
		powerTransmitted = powerTransmitted/(sD->xAnt*sD->yAnt);
		printf("z=%d\t energy=%f\n", z, powerTransmitted);
	}

	printf("propagation end \n");
}

/************************************************/
int main(int argc, char* argv[]) 
{
	sampData sD;
	modelData mD;

	/* Sample points for all model */	
	sD.xAnt = 1024; /* first transverse direction  */
	sD.yAnt = 1024; /* second transverse direction */
	sD.zAnt = 1024; /* depth                       */
	sD.dx = 0.1;    /* depth step                  */
	sD.dy = 0.1;
	sD.dz = 0.1;
		
	
	/* Size of box containing one bloodcell. Let this be 1/8 */
	/* of total number of samplepoints.                      */
	sD.xbox = sD->xAnt / 8;
	sD.ybox = sD->xAnt / 8;
	sD.zbox = sD->xAnt / 8;

	sD.zout = 128;
	sD.xout = 128;
	sD.iAnt=2;
		
		
	/* Field data */
	/* RBC max size 7.76 um is equivalent to sD.xbox */
	/* lambda = 632.8 nm = (sD.xbox/7.76)*0.6328 s.p */
	lambda = (sD.xbox/7.76)*0.6328; /* Sample points */
	mD.afreq = 2*Pi/lambda;          /* angular frequency c^{-1} s^{-1} */
	mD.epsilonRe = 1.977;            /* permittivity */
	mD.epsilonIm = 0.0002;           /* permittivity */
	mD.backRe = 1.809;
	mD.backIm = 0.0;

	/* Model contains 2D layer with ALL samplepoints */
	propagate( &sD, &mD );

	free(sD);
	free(mD);
}
/**********************************************************/