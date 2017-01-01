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

typedef struct bloodData{
	int amp, xc, yc, zc, dxx, dyy;
	double theta, fi, psi;
} bloodData;


typedef struct sampData{
	int xAnt, yAnt, zAnt, iAnt, xout, zout;
	double dx, dy, dz, fraction;
	int xbox, ybox, zbox;
} sampData;


/* Describes medium */
typedef struct modelData{
	double afreq; /* angular frequency */
	double epsilonRe, epsilonIm; /* inclusion permittivity */
	double backRe, backIm; /* backgound permittivity */
} modelData;


/* extern variables */
FILE *fphelp; /* debugging output file */
FILE *fphelp2; /* debugging output file */

/* output of structures ************************************/

void printfsampData(sampData *sD){
printf("Discretization data\n");
printf("xAnt = %d\t, yAnt = %d\t, zAnt = %d\n",sD->xAnt, sD->yAnt, sD->zAnt);
printf("dx = %.4f\t, dy = %.4f\t, dz = %.4f\n\n",sD->dx, sD->dy, sD->dz);
}


void fprintfsampData(FILE *fp, sampData *sD){
	fprintf(fp, "Discretization data\n");
	fprintf(fp, "xAnt = %d\t, yAnt = %d\t, zAnt = %d\n",
	sD->xAnt, sD->yAnt, sD->zAnt);
	fprintf(fp, "dx = %.4f\t, dy = %.4f\t, dz = %.4f\n\n",
	sD->dx, sD->dy, sD->dz);
}


void printfmodelData(modelData *mD){
	printf("Modelling data\n");
	printf("angular frequency = %.4f\n", mD->afreq);
	printf("inclusion = %.4f-i %.4f\n", mD->epsilonRe, mD->epsilonIm);
	printf("backgroud = %.4f-i %.4f\n", mD->backRe, mD->backIm);
}


void fprintfmodelData(FILE *fp, modelData *mD){
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

	/* Assert if out of boundaries */
	assert(sizeof(*d) == sizeof(double));
	assert(sizeof(*c) == sizeof(double));

}


void compmult(double ar, double ai, double br, double bi, double *cr,
double *ci)
{
	*cr = ar*br - ai*bi;
	*ci = ar*bi + ai*br;

    /* Assert if out of boundaries */
	assert(sizeof(*cr) == sizeof(double));
	assert(sizeof(*ci) == sizeof(double));
}


void compdiv(double ar, double ai, double br, double bi, double *cr, double
*ci)
{
	*ci = br*br + bi*bi;
	*cr = (ar*br + ai*bi)/(*ci);
	*ci = (br*ai - bi*ar)/(*ci);

	/* Assert if out of boundaries */
	assert(sizeof(*cr) == sizeof(double));
	assert(sizeof(*ci) == sizeof(double));
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
 * Funktionen allokerar upp en grid "dim1 x dim2 x dim3" med bloodData i
 * varje punkt.
 */
bloodData ****multialloc3blood(int dim1, int dim2, int dim3)
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
 * Skapar en 3D grid med samplingspunkter, varje lager slumpmässigt förskjutet
 * i förhållande till varandra.
 */
void initialize3DGrid(sampData *sD, bloodData ****bD, double wanted, FILE *fpmodel)
{
	int xmax=sD->xAnt, zmax=sD->zAnt, ymax=sD->yAnt;
	int xboxAnt, zboxAnt, yboxAnt;
	int zb, xb, yb;
	int xbox=sD->xbox, ybox=sD->ybox, zbox=sD->zbox; 
	
	
	xboxAnt = xmax/xbox;
	zboxAnt = zmax/zbox;
	yboxAnt = ymax/ybox;

	printf("Skapar skelett\n");
	/* Loops through all boxex, notice loops from 0->(xboxAnt-1) */
	for (xb=0; xb<xboxAnt; xb++){
		for (yb=0; yb<yboxAnt; yb++){
			for (zb=0; zb<zboxAnt; zb++) {
				bD[zb][yb][xb] -> zc = (zb+1)*zbox-zbox/2;
				bD[zb][yb][xb] -> xc = (xb+1)*xbox-xbox/2;
				bD[zb][yb][xb] -> yc = (yb+1)*ybox-ybox/2;

        /* Förskjut lager i x-y led */
				bD[zb][yb][xb] -> dxx = (int)random1(zbox);
				bD[zb][yb][xb] -> dyy = (int)random1(zbox);
				bD[zb][yb][xb] -> theta = random1(2*pi);
				bD[zb][yb][xb] -> fi = random1(2*pi);
				bD[zb][yb][xb] -> psi = random1(2*pi);
				bD[zb][yb][xb] -> amp = 1;

			}
		}
	}
	
	fpmodel = fopen( "model.dat", "w");
	fprintf(fpmodel, "%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n" , "x", "y", "z", "Amplitude", "Theta", "Fi","Psi", "xc", "yc","zc");
	for (zb=0; zb<zboxAnt; zb++){
		for (yb=0; yb<yboxAnt; yb++){
			for (xb=0; xb<xboxAnt; xb++) {
				fprintf(fpmodel,"%d\t %d\t %d\t %d\t %f\t %f\t %f\t %d\t %d\t %d\n", xb, yb, zb, bD[zb][yb][xb]->amp, bD[zb][yb][xb]->theta,
				bD[zb][yb][xb]->fi, bD[zb][yb][xb]->psi, bD[zb][yb][xb]->xc, bD[zb][yb][xb]->yc,
				bD[zb][yb][xb]->zc);
			}
		}
	}
}


/* 
* Fyll ett tvådimensionellt lager i modellen för fixt z!) 
*/

void fill2DLayer(sampData *sD, bloodData ***bD, MatInt& eps, int z){

	int xmax=sD->xAnt, zmax=sD->zAnt, ymax=sD->yAnt;
	int x, y, zb, xb, dxxx, dyyy,
	yb, xbox=sD->xbox, ybox=sD->ybox, zbox=sD->zbox, ix,
	jz, iy, antal=0;
	double a1, a2, a3, cost, sint, a=-0.00000005692,
	b=0.00000005692, abc, c=0.00227,
	sinf, cosf, sinp, cosp;

	for (y=0;y<ymax;y++) {
		for (x=0;x<xmax;x++) {

			/* Beräkna vilken "låda" vi är i */
			zb=z/zbox;
			xb=x/xbox;
			yb=y/ybox;

			/* Beräkna vinklar */
			cost = cos(bD[yb][xb]->theta);
			sint = sin(bD[yb][xb]->theta);
			sinf = sin(bD[yb][xb]->fi);
			cosf = cos(bD[yb][xb]->fi);
			cosp = cos(bD[yb][xb]->psi);
			sinp = sin(bD[yb][xb]->psi);

			ix = x - bD[yb][xb]->xc ;
			jz = z - bD[yb][xb]->zc ;
			iy = y - bD[yb][xb]->yc ;

			/* Transform coordinates */

			a1 = ix*(cosp*cosf-cost*sinf*sinp) + iy*(cosp*sinf+cost*cosf*sinp) + jz*sinp*sint;
			a2 = ix*(-sinp*cosf-cost*sinf*cosp) + iy*(-sinp*sinf+cost*cosf*cosp) + jz*cosp*sint;
			a3 = ix*(sint*sinf) + iy *(-sint*cosf) + jz*(cost);

			abc = a*(a1*a1+a2*a2)+b*(a1*a1+a2*a2)*(a1*a1+a2*a2)+c*a3*a3;
			if (abc < 1.0) {
				dxxx = bD[zb][1]->dxx;
				dyyy = bD[zb][1]->dyy;
				eps[(x+dxxx)%xmax][(y+dyyy)%ymax] = bD[yb][xb]-> amp;
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
	{ /*interpolation slowness */
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

	/* Obs: längd e nu "ymax x 2*xmax" */
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
void propagate(sampData *sD, modelData *mD, double wanted, FILE *fptrans, FILE *fpfield) {

	int xmax=sD->xAnt, zmax=sD->zAnt, ymax=sD->yAnt, zout=sD->zout, xout=sD->xout, imax=2;
	double dx=sD->dx, dz=sD->dz, dy=sD->dy;
	int xbox=sD->xbox, ybox=sD->ybox, zbox=sD->zbox;
	int x, y, z;
	double energytrans=0.0, epsr=mD->epsilonRe;
	extern FILE *fphelp, *fphelp2;
    
    int a = 0;
    MatDoub umig(ymax, 2*xmax, a);
    MatDoub uref(ymax, 2*xmax, a);
    
    /* Allocate memory for layers */
    MatInt model(ymax, xmax, a);
    MatInt nextmodel(ymax, xmax, a);

	bloodData ****bD;

	/* Allocate memory for 3D grid, every point containing bloodData.
	 * This 3D grid contains ONLY centers for bloodcells together with
	 * their properties.
	 */

	bD = multialloc3blood(zmax/zbox, ymax/ybox, xmax/xbox);
	wanted=0.4;

	/* Initialize the 3D grid, blood cell centers only */
	initialize3DGrid(sD, bD, wanted, fphelp);

	/* Set boundary field with initiale values */
	for (y=0; y<ymax; y++) {
		for (x=0; x<2*xmax; x+=2) {
			umig[y][x] = 1000; //1.0;
			umig[y][x+1] = 0.0;
			
			/* Reflected field */
			uref[y][x] = 0;
			uref[y][x+1] = 0;
		}
	}

	printf("Propagate \n");
	fphelp = fopen("powerout", "w");

	/* Start propagate from z=0 to z=zmax-1 */
	for (z=0; z<zmax; z++) {
		
		/* Fill two layers in xy-plane with sample points */ 
		fill2DLayer(sD, bD[z/zbox], model ,z);
		fill2DLayer(sD, bD[(z+1)/zbox], nextmodel,z+1);
        
        printf("Kalle\n");

		/* Migrate between the two layers */
		migrate(sD, mD, model, umig);
        
		/* Calculate interaction between two layers */
		interaction(sD, mD, model, nextmodel, umig, uref);

		/* Calculate power in z by summing over x-y plane */
		for (y=0; y<ymax; y++) {
			for (x=0; x<2*xmax; x+=2) {
				energytrans += (umig[y][x]*umig[y][x]+umig[y][x+1]*umig[y][x+1])*(1.0+epsr*model[y][x/2]);
			}
		}

			
		/* normalize */
				energytrans = energytrans/(xmax*ymax);
				fprintf(fphelp,"%f\n", energytrans);
				printf("z=%d\t energy=%f\n", z, energytrans);
				energytrans = 0.0;


		//if (fptrans) {
		//	for (y=0; y<ymax; y++) {
		//		for (x=0; x<2*xmax; x++){
		//			fprintf(fptrans, "%lf\n", umig[y][x]);
		//		}
		//	}
		//}
	}


	fclose(fptrans);
	fclose(fphelp);
	printf("i7 %d\n");
	/* multidealloc3(v1, imax, ymax, 2*xmax);
	multidealloc2(umig, ymax, 2*xmax); */
	//free(uref);
	//multidealloc2int(nextmodel, sD->yAnt, sD->xAnt);
	printf("propagation end \n");
}
/* the main: input of data....
************************************************/
int main(int argc, char* argv[]) {
	int x, y, z, r, error=0, error1=0;
	double tid, wanted=0.4;
	char
	outfile[80], /* output name */
	logfile[80], /* program log file */
	transfile[80], /* transmitted field file */
	fieldfile[80], /* internal field file */
	helpfile[80], helpfile2[80]; /* help files */
	FILE *fptrans, *fpfield, *fplog; /* output files */

	extern FILE *fphelp, *fphelp2;
	time_t starttid = time(NULL), sluttid;
	int t1, t2;
	sampData *sD;
	modelData *mD;
	sD = (sampData *) malloc(sizeof(*sD));
	mD = (modelData *) malloc(sizeof(*mD));


	/* Sample points for all modell */	
	sD->xAnt = 1024;      /* first transverse direction */
	sD->yAnt = 1024; /* second transverse direction */
	sD->zAnt = 1024; //256;      /* depth */
	sD->dx = 0.1;          
	sD->dy = 0.1; 
	sD->dz = 0.1;          /* depth step */
		
	
	/* Size of box containing one bloodcell.     */
	/* Let this be 1/8 of total size.            */
	/* Having 128x128 gridpoints then implies    */
	/* 8x8 boxes each containing 16^3 = 4096     */
	/* samplepoints.                             */
	sD->xbox = sD->xAnt / 8;
	sD->ybox = sD->xAnt / 8;
	sD->zbox = sD->xAnt / 8;

	sD->zout = 128; //atoi(argv[11]);
	sD->xout = 128;
	sD->iAnt=2; /* pre defined quanitites */
		
		
	/* Field data */
	mD->afreq = 0.1; //atof(argv[8]); /* frequency */
	mD->epsilonRe = 1.94; //atof(argv[9]); /*permittivity */
	mD->epsilonIm = 0.0002; //atof(argv[10]); /*permittivity */
	mD->backRe = 1.81;
	mD->backIm = 0.0;
		

	strcpy( logfile, outfile);
	strcat( logfile, ".log");
	fplog = fopen( logfile, "w");
	strcpy( transfile, outfile);
	strcat( transfile, "_trans.out");
	fptrans = fopen( transfile, "w");
	strcpy( fieldfile, outfile);
	strcat( fieldfile, "_field.out");
	fpfield = fopen( fieldfile, "w");
	strcpy( helpfile, outfile);
	strcat( helpfile, "_help.out");
	fphelp = fopen( helpfile, "w");
	strcpy( helpfile2, outfile);
	strcat( helpfile2, "_help2.out");
	fphelp2 = fopen( helpfile2, "w");
	printfsampData(sD);
	printfmodelData(mD);

	t1 = clock();
	/* input of the model */
	x = 0;
	z = 0;
	r = 0;
	y = 0;
	error1 = 0;
	/* Model contains 2D layer with ALL samplepoints */
	if (error==0) {
		printf("2 propagate\n");
		propagate(sD, mD, wanted, fptrans, fpfield);
	} 
	else
	printf("Terminating program, Input Error: %d.\n", error);
	sluttid = time(NULL);
	fprintf(fplog, "end on ");
	fprintf(fplog, ctime(&sluttid));
	tid = difftime(sluttid, starttid);
	fprintf(fplog, "total time: %.0f\t: %.0f:%.0f\n", tid, tid/3600, tid/60);
	t2 = clock();
	fprintf(fplog, "CPU clock: %d\t%d\t%e\n", t1, t2, (t2-t1)/1000000.0);

	free(sD);
	free(mD);
}
/**********************************************************/