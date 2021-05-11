/* Bo-Göran Wallner */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <random>
#include <complex>

//#include <time.h>
#include <ctime>
#include <assert.h>

#include "nr3.h"
#include "fourier_ndim.h"

/* Some constants*/
#define twopi                  6.283185307177959
#define pi                     3.14159265358
#define max32bit               4294967295
#define initFieldValue         1000.0

/* Development environment */
#define GNU_LINUX              0

/* Transversal size 2^N */
#define MODEL_TRANSVERSAL_SIZE 1024
#define MODEL_DEPTH_SIZE       1024

/* Choose one angle setup for different simulations */
#define ANGLES_RANDOM          0
#define ANGLE_THETA_ZERO       0
#define ANGLE_THETA_PI_HALF    0
#define ANGLE_THETA_PI_FOURTH  1

/* Enable backward field */
#define ENABLE_BREMMER_REFLECTION 1

/* Number of layers to simulate in propagation direction */
#define SIMULATION_DEPTH       1020

/**** Flags for debug output in textfiles ****/

/* Print multipe layer of RBC to file */
#define MULTIPLE_LAYER_NUMBER_TO_FILE       128
/* Where to start writing geometry */
#define START_MULTIPLE_LAYER_NUMBER_TO_FILE 510
/* Print all BloodData to file */
#define BLOOD_DATA_TO_FILE                  0
/* Print transmitted intensity */
#define TRANS_INTENSITY_TO_FILE             1
/* Print vars in propagate() */
#define PRINT_MIGRATE_DATA                  0           
/* Print field absolute to file */
#define PRINT_MIGRATED_FIELD_TO_FILE        0
/* Print field every every N:th layer */
#define PRINT_FIELD_EVERY_NTH_LAYER         0

/* File handles */
static FILE *fpSampleLayer;
static FILE *fpBloodData;
static FILE* fpTransPower;
static FILE* fpInteractionData;
static FILE* fpMigrateData;
static FILE* fpMigratedField;

/* RBC width in sample points */
static unsigned int rbcWidthInSampPoints;

/* Counters for determing total number of    */
/* samplepoints within RBC and in Background */
static unsigned long nbrOfSamplespointsInRbc;
static unsigned long nbrOfSamplespointsInBackground;

typedef struct
{
    int amp, xc, yc, zc, dx, dy;
    double theta, fi, psi;
} bloodData;

typedef struct 
{
    int xAnt, yAnt, zAnt, iAnt;
    double dx, dy, dz;
    int xbox, ybox, zbox;
} sampData;


/* Describes medium */
typedef struct
{
    double afreq;                /* angular frequency */
    double epsilonRe, epsilonIm; /* inclusion permittivity */
    double backRe, backIm;       /* background permittivity */
} modelData;


/* output of structures ************************************/

/* Function prints the transmitted intensity */
static void fprintftransIntensityData(double transPower)
{
    char buf[40];
    snprintf(buf, sizeof(buf), "data/fieldData/transpowerflux.txt");

#if ( GNU_LINUX == 1 )
    fpTransPower = fopen(buf, "a");
#else
    fopen_s(&fpTransPower, buf, "a");
#endif

    fprintf(fpTransPower, "Transmitted intensity = %.8f\n", transPower);
    fclose(fpTransPower);
}


/* Function prints one geometry layer to file */
static void fprintfsampData(sampData *sD, MatInt& geometry2D, int z)
{
    int x, y;
    char buf[40];
    snprintf(buf, sizeof(buf), "data/geometryData/samplelayer%d.txt", z);

#if ( GNU_LINUX == 1 )
    fpSampleLayer = fopen(buf, "w");
#else
    fopen_s(&fpSampleLayer, buf, "w");
#endif

    for (y = 0; y < sD->yAnt; y++)
    {
        for (x = 0; x < sD->xAnt; x++)
        {
            fprintf(fpSampleLayer, "%d\t", geometry2D[y][x]);
        }
        fprintf(fpSampleLayer, "\n");
    }
    fclose(fpSampleLayer);
}

/* Function prints one migrated field data to file */
static void fprintfMigratedFieldData(sampData* sD, MatDoub& field, int z)
{
    int x, y;
    double absoluteField;
    char buf[40];
    snprintf(buf, sizeof(buf), "data/fieldData/migratedfield%d.txt", z);

#if ( GNU_LINUX == 1 )
    fpMigratedField = fopen(buf, "w");
#else
    fopen_s(&fpMigratedField, buf, "w");
#endif

    for (y = 0; y < sD->yAnt; y++)
    {
        for (x = 0; x < sD->xAnt; x++)
        {
            absoluteField = field[y][2*x]*field[y][2*x] + field[y][2*x+1]*field[y][2*x+1];
            fprintf(fpMigratedField, "%.8f\t", sqrt(absoluteField)/ initFieldValue);
        }
        fprintf(fpMigratedField, "\n");
    }
    fclose(fpMigratedField);
}

static void fprintMigrateData(double epsr, double epsi, double w, double eta, double kr, 
                              double ki, double k2r, double ga2i)
{
    char buf[40];
    snprintf(buf, sizeof(buf), "data/modelData/migratedata.txt");

#if ( GNU_LINUX == 1 )
    fpMigrateData = fopen(buf, "a");
#else
    fopen_s(&fpMigrateData, buf, "a");
#endif

    fprintf(fpMigrateData, "permitivity = %.8f-i %.15f\n", epsr, epsi);
    fprintf(fpMigrateData, "angular freq = %.8f\t eta = %.8f\t kr = %.8f\t ki = %.8f\n", w, eta, kr, ki);
    fprintf(fpMigrateData, "k2r = %.8f\t ga2i = %.15f\n", k2r, ga2i);
    fprintf(fpMigrateData, "\n");
    fclose(fpMigrateData);
}

void fprintfmodelData(FILE *fp, modelData *mD)
{
    fprintf(fp, "Modelling data\n");
    fprintf(fp, "angular frequency = %.4f\n", mD->afreq);
    fprintf(fp, "inclusion = %.4f-i %.4f\n", mD->epsilonRe, mD->epsilonIm);
    fprintf(fp, "background = %.4f-i %.4f\n", mD->backRe, mD->backIm);
}

static void fprintfBloodData( bloodData**** bD, int zdim, int ydim, int xdim )
{
    char buf[40];
    int x, y, z;
    snprintf(buf, sizeof(buf), "data/modelData/bloodData.txt");

#if ( GNU_LINUX == 1 )
    fpBloodData = fopen(buf, "w");
#else
    fopen_s(&fpBloodData, buf, "w");
#endif

    fprintf(fpBloodData, "**** BloodData 3D Geometry ****\n");
    fprintf(fpBloodData, "\n");
    fprintf(fpBloodData, "\n");
    fprintf(fpBloodData, "xc = center x coordinate\n");
    fprintf(fpBloodData, "yc = center y coordinate\n");
    fprintf(fpBloodData, "zc = center z coordinate\n");
    fprintf(fpBloodData, "th = Euler theta angle (deg)\n");
    fprintf(fpBloodData, "fi = Euler fi angle (deg)\n");
    fprintf(fpBloodData, "psi = Euler psi angle (deg)\n");
    fprintf(fpBloodData, "ry = random y displacement\n");
    fprintf(fpBloodData, "rx = random x displacement\n");

    fprintf(fpBloodData, "\n");
    fprintf(fpBloodData, "\n");

    for (z = 0; z < zdim; z++)
    {
        for (y = 0; y < ydim; y++)
        {
            for (x = 0; x < xdim; x++)
            {
                fprintf(fpBloodData, "zc=%d\t", bD[z][y][x]->zc);
                fprintf(fpBloodData, "yc=%d\t", bD[z][y][x]->yc);
                fprintf(fpBloodData, "xc=%d\t", bD[z][y][x]->xc);
                fprintf(fpBloodData, "th=%d\t", (int)(bD[z][y][x]->theta * 180.0/pi));
                fprintf(fpBloodData, "fi=%d\t", (int)(bD[z][y][x]->fi * 180.0/pi));
                fprintf(fpBloodData, "psi=%d\t", (int)(bD[z][y][x]->psi * 180.0/pi));
                fprintf(fpBloodData, "ry=%d\t", bD[z][y][x]->dy);
                fprintf(fpBloodData, "rx=%d\t", bD[z][y][x]->dx);
                fprintf(fpBloodData, "\n");
            }
            fprintf(fpBloodData, "\n");
        }
        fprintf(fpBloodData, "\n");
    }
    fclose(fpBloodData);
}

/* some complex valued algebra *****************************/
static void compsqrt(double a, double b, double *c, double *d)
{
    complex<double> dval;
    dval = std::sqrt(std::complex<double>(a, b));
    *c = dval.real();
    *d = dval.imag();
}

/* Multiply (ar + i*ai)(br + i*bi) */
static void compmult(double ar, double ai, double br, double bi, double *cr, double *ci)
{
    complex<double> dval;
    dval = (std::complex<double>(ar, ai)) * (std::complex<double>(br, bi));
    *cr = dval.real();
    *ci = dval.imag();
}

static void compdiv(double ar, double ai, double br, double bi, double *cr, double *ci)
{
    complex<double> dval;
    dval = (std::complex<double>(ar, ai)) / (std::complex<double>(br, bi));
    *cr = dval.real();
    *ci = dval.imag();
}

static double compabs(double ar, double ai)
{
    return std::abs(std::complex<double>(ar, ai));
}

static double absolut(double a)
{
    return (a < 0) ? -a : a;
}

static double random1(double a)
{
    std::random_device dev;
    std::mt19937 engine(dev());
    std::uniform_real_distribution<double> dist{ 0.0, 1.0 };
    return (a * dist(engine));
}

/* Calculate FFT */
static void fourn_wrapper(MatDoub& u1, int xmax, int ymax, bool inverse)
{
    int i, y, x;
    bool stop;

    /* Create a Vector with double and all zeros */
    VecDoub v1(2 * xmax * ymax, 0.0);

    /* 2D input to NR fourn() */
    VecInt nn(3);
    nn[0] = 1;
    nn[1] = ymax;
    nn[2] = xmax;


    /* Copy u1[ymax][2*xmax] -> v1[2*xmax*ymax] */
    i = 0;
    for (y = 0; y < ymax; y++)
    {
        for (x = 0; x < (2*xmax-1); x += 2)
        {
            v1[i] = u1[y][x];
            v1[i + 1] = u1[y][x + 1];
            i += 2;

            /* DEBUG test */
            if (v1[i] > initFieldValue || v1[i + 1] > initFieldValue)
            {
                stop = true;
            }
        }
    }

    if (false == inverse)
    {
        /* Call FFT */
        fourn(v1, nn, 1);
    }
    else
    {
        /* Inverse FFT */
        fourn(v1, nn, -1);
    }

    /* Copy back contents from v1[2*xmax*ymax] -> u1[ymax][2*xmax] */
    i = 0;
    for (y = 0; y < ymax; y++)
    {
        for (x = 0; x < (2*xmax-1); x += 2)
        {
            if( true == inverse )
            {
                /* Inverse transform needs division with all dimensions */
                u1[y][x] = v1[i] / (xmax * ymax);
                u1[y][x + 1] = v1[i + 1] / (xmax * ymax);

                /* DEBUG test */
                if (u1[y][x] > initFieldValue || u1[y][x+1] > initFieldValue)
                {
                    stop = true;
                }
            }
            else
            {
                u1[y][x] = v1[i];
                u1[y][x + 1] = v1[i + 1];
            }
            i += 2;
        }
    }
}

/* 
 * Memoryallocation of bloodData 3D array
 */
static bloodData ****mallocBloodData3DArray(int dim1, int dim2, int dim3)
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
static void initBloodData3DArray(sampData *sD, bloodData ****bD)
{
    int nbrOfXbox, nbrOfYbox, nbrOfZbox;
    int zb, xb, yb, dx, dy;

    /* Calculate number of cells. E.g if we use 1024 sample points    */
    /* and each cell is 128 points we get 1024/128 = 8 cells          */
    nbrOfXbox = sD->xAnt/sD->xbox;
    nbrOfZbox = sD->zAnt/sD->zbox;
    nbrOfYbox = sD->yAnt/sD->ybox;

    for (zb = 0; zb < nbrOfZbox; zb++)
    {
        /* Displace centers of whole xy-plane layer         */
        /* in y-direction to create randomness in geometry. */
        dy = (int)random1((double)sD->ybox);

        for (yb = 0; yb < nbrOfYbox; yb++)
        {
            /* For every fix yb translate in x-direction */
            dx = (int)random1((double)sD->xbox);
            for (xb = 0; xb < nbrOfXbox; xb++)
            {
                bD[zb][yb][xb]->dx = dx;
                bD[zb][yb][xb]->dy = dy;
                bD[zb][yb][xb]->zc = (zb+0.5)*sD->zbox;
                bD[zb][yb][xb]->xc = (xb+0.5)*sD->xbox;
                bD[zb][yb][xb]->yc = (yb+0.5)*sD->ybox;
                bD[zb][yb][xb]->amp = 1;

                /* E.g. if we have 1024 samplepoints and choose 128 size cells        */
                /* we get 8 cells having centers (0+0.5)*128 = 64, (1+0.5)*128 = 192, */
                /* (2+0.5)* 128 = 320, ... (7+0.5)*128 = 960                          */

                /* Assign values for the Euler angles */
#if    (ANGLES_RANDOM == 1)
                bD[zb][yb][xb]->theta = random1(pi);
                bD[zb][yb][xb]->fi = random1(2*pi);
                bD[zb][yb][xb]->psi = random1(2*pi);
#endif

#if    (ANGLE_THETA_ZERO == 1)
                bD[zb][yb][xb]->theta = 0;
                bD[zb][yb][xb]->fi = random1(2*pi);
                bD[zb][yb][xb]->psi = random1(2*pi);
#endif

#if    (ANGLE_THETA_PI_HALF == 1)
                bD[zb][yb][xb]->theta = pi/2;
                bD[zb][yb][xb]->fi = random1(2*pi);
                bD[zb][yb][xb]->psi = random1(2*pi);
#endif

#if    (ANGLE_THETA_PI_FOURTH == 1)
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
            /* See e.g. Goldstein page 146-147 x'=A*x=BCD*x (eq 4-46)   */
            xe = xl*(cosp*cosf-cost*sinf*sinp) + yl*(cosp*sinf+cost*cosf*sinp) + zl*(sinp*sint);
            ye = xl*(-sinp*cosf-cost*sinf*cosp) + yl*(-sinp*sinf+cost*cosf*cosp) + zl*(cosp*sint);
            ze = xl*(sint*sinf) + yl *(-sint*cosf) + zl*(cost);

            /* So far we have not assumed anything about the exact     */
            /* geometry. We have barely stated, we have any object     */
            /* shape in (x,y,z) and translate this object into a       */
            /* local coordinate system (xe,ye,ze) by translation       */
            /* and Euler rotation (since we decided to describe the    */
            /* object with Euler angles).                              */

            /******************* BICONCAVE DISK / BLOODCELL APPROACH ***********************/

            /* Suggestion by A.Karlsson "Numerical Simulations of Light Scattering          */
            /* by Red Blood Cells" the equation for the biconcave disk can be               */
            /* described as: D(x) = Sqrt(1-(x/R0)^2) * (C0+C2*(x/R0)^2+C4*(x/R0)^4))        */
            /* where C0=0.81um, C2=7.83um, C4=-4.39um and R0=3.91um and this corresponds    */
            /* to a biconcave disk with diameter a=7.76um and max thickness b=2.55um.       */
            /* Thus, the a=7.76um shall correspond to rbcWidthInSampPoints samplepoints.    */
            /* We apply the following tranformations: x -> Sqrt(xe^2+ye^2) and for D(x)     */
            /* we set D(x)^2 -> (2ze)^2 = 4*ze^2. This give a 3D expression and we can use  */
            /* 4*ze^2 - (1-(Re/R0)^2) * (C0+C2*(Re/R0)^2+C4*(Re/R0)^4))^2 > 0 and if        */
            /* this expression holds then the ze is *not* within the biconcave disk.        */
            /* Also Re=sqrt(xe^2+ye^2) < half cell size to determine if (xe,ye) is within   */
            /* the disk or not.                                                             */

            /* NOTICE: One could add more decimals to C0, C2 and C4 to get a more exact     */
            /* condition. Now the boundary can differ slightly for the RBC making it to     */
            /* take slightly more space than the size of the "box". In practice this is a   */
            /* neglectible problem.                                                         */

            R0 = 64.4948;
            C0 = 13.3608;
            C2 = 129.1546;
            C4 = -72.41243;
            /* Calculate (xe,ye) corresponding radius */
            Re = sqrt(pow(xe,2) + pow(ye,2));

            /* The whole xy-layer shall be translated */
            dx = bD[yb][xb]->dx;
            dy = bD[yb][xb]->dy;

            /* Check if (xe,ye) is within disk since expression only valid within */
            /* it and that exclude e.g. corners of the cell.                      */
            if(Re < rbcWidthInSampPoints/2)
            {
                /* Re is valid to use in the expression */
                D_xDiff = 4*pow(ze,2) - (1-pow(Re/R0,2)) * pow(C0 + C2*pow(Re/R0,2) + C4*pow(Re/R0,4), 2);
                /* Check if ze is outside or within the biconcave disk */
                if ( D_xDiff > 0.0)
                {
                    /* ze is outside the disk */
                    geometry2D[(x+dx)%sD->xAnt][(y+dy)%sD->yAnt] = 0;
                    nbrOfSamplespointsInBackground++;
                }
                else
                {
                /* ze is within the biconcave disk */
                    geometry2D[(x+dx)%sD->xAnt][(y+dy)%sD->yAnt] = 1;
                    nbrOfSamplespointsInRbc++;
                }
            }
            else
            {
               /* Re is outside the disk */
                geometry2D[(x+dx)%sD->xAnt][(y+dy)%sD->yAnt] = 0;
                nbrOfSamplespointsInBackground++;
            }
        }
    }

#if ( MULTIPLE_LAYER_NUMBER_TO_FILE > 0)
    /* Print layers with start START_MULTIPLE_LAYER_NUMBER_TO_FILE to files */
    if ((z > START_MULTIPLE_LAYER_NUMBER_TO_FILE) && 
        (z < (MULTIPLE_LAYER_NUMBER_TO_FILE + START_MULTIPLE_LAYER_NUMBER_TO_FILE)))
    {
        fprintfsampData(sD, geometry2D, z);
    }
#endif

}

/* migrate the fields one layer ****************************************/
static void migrate(sampData *sD, modelData *mD, MatInt& slow, MatDoub& u1)
{
    int xmax=sD->xAnt, zmax=sD->zAnt, imax=sD->iAnt, ymax=sD->yAnt;
    double dx=sD->dx, dz=sD->dz, dy=sD->dy ;
    double w=mD->afreq, eta;
    double s, xiX, xiY, smin, smax, b, temp, propr, propi;
    double ii;
    double epsr, epsi, ga2r, ga2i, gar, gai, kr, ki, k2r;
    int x, i, y, pr=0, ev=0;
    int imMax, reMax;

    /* Create sD.iAnt number of vectors */
    vector<MatDoub> v1(sD->iAnt);
    for( i=0; i<sD->iAnt; i++ )
    {
        /* Resize each vector to contain the field */
        v1[i].resize(ymax, 2*xmax);
        v1[i].assign(ymax, 2 * xmax, 0.0);
    }

    /* Smooth out the 2D geometry */
    smin=1000.0;
    smax=0.0;
    for ( y=0; y<ymax; y++ ) 
    {
        for ( x=0; x<xmax; x++ )
        { /* slowness max and min */
            temp = slow[y][x];
            smax = (smax > temp) ? smax : temp;
            smin = (smin < temp) ? smin : temp;
        }
    }

    /* Calculate FFT for the input field */
    fourn_wrapper( u1, xmax, ymax, false );

    eta=0.0; /* is OK with 0.0 */
    for ( i=0; i<imax; i++ )
    {   /*interpolation slowness */
        s = smin + i*(smax-smin)/(imax-1);                                    /* slowness */
        epsr = (mD->backRe + s * (mD->epsilonRe - mD->backRe)) / mD->backRe;  /* Re permittivity, min = mD->back, max = mD->epsilonRe */
        epsi = -s* mD->epsilonIm;                                              /* Im permittivity*/
        compmult(eta, w, eta, w, &kr, &ki);                                   /* wave number, c^{-1} s */
        compmult(kr, ki, epsr, epsi, &k2r, &ga2i);                            /* sqr wave number, (kr+i*ki)(eps + i*epsi), .^{2} */

#if ( PRINT_MIGRATE_DATA == 1 )
        fprintMigrateData(epsr, epsi, w, eta, kr, ki, k2r, ga2i);
#endif

        /* For ordering from fourn() see Numerical recipes users guide   */
        /* which gives separate cases of y=0 and y=ymax/2. All comments  */
        /* given for xmax=ymax=1024 but calculations not restricted to   */
        /* these ranges.                                                 */

        /* Start indexes for array-indexes going "backwards" */
        imMax = 2*xmax-1; /* If xmax=1024 -> 2047 */
        reMax = 2*xmax-2; /* If xmax=1024 -> 2046*/

        /************* Handle the 0 frequency ************/
        y = 0;
        xiY = y/(ymax*dy)*twopi;       /* transverse wavenumber */
        /* x=0..1022 */
        for (x = 0; x < xmax; x += 2)
        {
            xiX = x / (xmax * 2.0 * dx) * twopi; /* transverse wavenumber */
            ga2r = k2r + xiX * xiX + xiY * xiY;
            compsqrt(ga2r, ga2i, &gar, &gai);
            propr = cos(-dz * gai) * exp(-dz * gar);
            propi = sin(-dz * gai) * exp(-dz * gar);
            /* xmax=1024 -> u[0][0], u[0][1] ... u[0][1022], u[0][1023] */
            compmult(propr, propi, u1[y][x], u1[y][x+1], &v1[i][y][x], &v1[i][y][x+1]);
            /* xmax=1024 -> u[0][2046], u[0][2047] ... u[0][1024], u[0][1025] */
            compmult(propr, propi, u1[y][reMax-x], u1[y][imMax-x], &v1[i][y][reMax-x], &v1[i][y][imMax-x]);
        }
        /********** END - Handle the 0 frequency *********/

        /********** Handle the mid +/- frequency *********/
        y = ymax/2;
        xiY = y/(ymax*dy)*twopi;       /* transverse wavenumber */
        /* x=0..1022 */
        for (x = 0; x < xmax; x += 2)
        {
            xiX = x / (xmax * 2.0 * dx) * twopi; /* transverse wavenumber */
            ga2r = k2r + xiX * xiX + xiY * xiY;
            compsqrt(ga2r, ga2i, &gar, &gai);
            propr = cos(-dz * gai) * exp(-dz * gar);
            propi = sin(-dz * gai) * exp(-dz * gar);
            /* xmax=1024 -> u[512][0], u[512][1] ... u[512][1022], u[512][1023] */
            compmult(propr, propi, u1[y][x], u1[y][x + 1], &v1[i][y][x], &v1[i][y][x + 1]);
            /* xmax=1024 -> u[512][2046], u[512][2047] ... u[512][1024], u[512][1025] */
            compmult(propr, propi, u1[y][reMax - x], u1[y][imMax - x], &v1[i][y][reMax - x], &v1[i][y][imMax - x]);
        }
        /******* END - Handle the mid +/- frequency ******/

        /******** Handle all other frequencys ************/
        /* y = 1..511 */
        for ( y=1; y<ymax/2; y++ )
        {
            xiY = y/(ymax*dy)*twopi;                 /* transverse wavenumber */

            /* x=0..1022 */
            for (x = 0; x < xmax; x += 2)
            {
                xiX = x / (xmax * 2.0 * dx) * twopi; /* transverse wavenumber */
                ga2r = k2r + xiX * xiX + xiY * xiY;
                compsqrt(ga2r, ga2i, &gar, &gai);
                propr = cos(-dz * gai) * exp(-dz * gar);
                propi = sin(-dz * gai) * exp(-dz * gar);

                /* Handle y array-index 1-511 */
                /* xmax=1024 -> u[y][0], u[y][1] ... u[y][1022], u[y][1023] */
                compmult(propr, propi, u1[y][x], u1[y][x+1], &v1[i][y][x], &v1[i][y][x+1]);
                /* xmax=1024 -> u[y][2046], u[y][2047] ... u[y][1024], u[y][1025] */
                compmult(propr, propi, u1[y][reMax-x], u1[y][imMax-x], &v1[i][y][reMax-x], &v1[i][y][imMax-x]);

                /* Handle y array-index ymax-1=1023...513. ymax/2=512 separately handled above! */
                /* xmax=1024 -> u[ymax-y][0], u[ymax-y][1] ... u[ymax-y][1022], u[ymax-y][1023] */
                compmult(propr, propi, u1[ymax-y][x], u1[ymax-y][x+1], &v1[i][ymax-y][x], &v1[i][ymax-y][x+1]);
                /* xmax=1024 -> u[ymax-y][2046], u[ymax-y][2047] ... u[ymax-y][1024], u[ymax-y][1025] */
                compmult(propr, propi, u1[ymax-y][reMax-x], u1[ymax-y][imMax-x], &v1[i][ymax-y][reMax-x], &v1[i][ymax-y][imMax-x]);
            }
        }

        /******** END - Handle all other frequencys ******/

        /* Make inverse FFT */
        fourn_wrapper(v1[i], xmax, ymax, true);

    } /* interpolation slowness */

    if ( ((smax-smin)>0.001) && (imax>1) )
    {
        for ( y=0; y<ymax; y++ ) 
        {
            /* xmax=1024 -> x=0,1,...,1023 -> max{2x+1} = 2*1023+1 = 2047  */
            for ( x=0; x<xmax; x++ )
            { 
                /* interpolation */
                b = modf((slow[y][x]-smin)/(smax-smin)*(imax-1), &ii);
                i = (int)ii;
                if ( i<(imax-1) )
                {
                    u1[y][2*x] = ((1.0-b)*v1[i][y][2*x] + b*v1[i+1][y][2*x]);
                    u1[y][2*x+1] = ((1.0-b)*v1[i][y][2*x+1] + b*v1[i+1][y][2*x+1]);
                }
                else 
                {
                    u1[y][2*x] = v1[i][y][2*x];
                    u1[y][2*x+1] = v1[i][y][2*x+1];
                }
            }
        }
    }
    else
    {
        for ( y=0; y<ymax; y++ ) 
        {
            /* Plane iteration over geometry. If xmax=1024 then x=0,1,...,2047. */
            for ( x=0; x<2*xmax; x++ )
            {  
                /*interpolation */
                u1[y][x] = v1[0][y][x];
            }
        }
    }
}

static void interaction( sampData *sD, modelData *mD,
                         MatInt &slow1, MatInt &slow2,
                         MatDoub& uin, MatDoub& urefl) {

    int xmax=sD->xAnt, zmax=sD->zAnt, ymax=sD->yAnt;
    double dx=sD->dx, dz=sD->dz, dy=sD->dy;
    double w=mD->afreq, eta=0.0;
    double epsr, epsi, kr, ki, k2r, gar, gai, ga2r, ga2i;
    double bsr, bsi, bkr, bki, bk2r, bgar, bgai, bga2r, bga2i;
    double refr, refi, xiX, xiY;
    int x, y, refl=0;
    int reMax, imMax;

    MatInt refp(ymax, xmax);
    /* Checks whether RBC->Background or Background->RBC */
    /* Essentially refl=1 most of the time.              */
    for (y=0; y<ymax; y++)
    {
        for (x=0; x<xmax; x++)
        {
            refp[y][x] = 0;
            if (slow1[y][x]-slow2[y][x] > 0.5)
            {
                /* RBC -> Background */
                refp[y][x] = 1;
                refl = 1;
            }
            else if (slow2[y][x]-slow1[y][x] > 0.5)
            {
                /* Background->RBC */
                refp[y][x] = -1;
                refl = 1;
            }
        }
    }

    /* Length is "ymax x 2*xmax" */
    if (refl == 1) 
    {
        /* Transform in-field */
        fourn_wrapper( uin, xmax, ymax, false );     /* spatial FFT */
        
        /* Set permitivity */
        epsr = (mD->backRe + 1.0 * (mD->epsilonRe - mD->backRe)) / mD->backRe; /* Re permittivity, min = mD->back, max = mD->epsilonRe */
        epsi = -mD->epsilonIm;                                                 /* Im permittivity*/

        /* Compute wavenumber for cells */
        compmult(eta, w, eta, w, &kr, &ki);           /* wavenumber, c^{-1} s */
        compmult(kr, ki, epsr, epsi, &k2r, &ga2i);    /* sqr wavenumber, .^{2} */
        
        /* Include background */
        bsr = 1.0;                                    /* Re of the background */
        bsi = 0.0;

        /* Compute wavenumber for background */
        compmult(bsr, bsi, eta, w, &bkr, &bki);       /* wavenumber, c^{-1} s */
        compmult(bkr, bki, bkr, bki, &bk2r, &bga2i);  /* sqr wavenumber, .^{2} */

        /* For ordering from fourn() see Numerical recipes users guide   */
        /* which gives separate cases of y=0 and y=ymax/2. All comments  */
        /* given for xmax=ymax=1024 but calculations not restricted to   */
        /* these ranges.                                                 */

        /* Start indexes for array-indexes going "backwards" */
        imMax = 2*xmax-1; /* If xmax=1024 -> 2047 */
        reMax = 2*xmax-2; /* If xmax=1024 -> 2046*/

        /************* Handle the 0 frequency ************/
        y = 0;
        xiY = y/(ymax*dy)*twopi;       /* transverse wavenumber */
        /* x=0..1022 */
        for (x = 0; x < xmax; x += 2)
        {
            xiX = x / (xmax * 2.0 * dx) * twopi; /* transverse wavenumber */
            ga2r = k2r + xiX * xiX + xiY * xiY;
            compsqrt(ga2r, ga2i, &gar, &gai);

            bga2r = bk2r + xiX*xiX + xiY*xiY;
            compsqrt(bga2r, bga2i, &bgar, &bgai);
            compdiv(gar-bgar, gai-bgai, gar+bgar, gai+bgai, &refr, &refi);

            /* xmax=1024 -> u[0][0], u[0][1] ... u[0][1022], u[0][1023] */
            compmult(refr, refi, uin[y][x], uin[y][x+1], &urefl[y][x], &urefl[y][x+1]);
            /* xmax=1024 -> u[0][2046], u[0][2047] ... u[0][1024], u[0][1025] */
            compmult(refr, refi, uin[y][reMax-x], uin[y][imMax-x], &urefl[y][reMax-x], &urefl[y][imMax-x]);
        }
        /********** END - Handle the 0 frequency *********/

        /********** Handle the mid +/- frequency *********/
        y = ymax/2;
        xiY = y/(ymax*dy)*twopi;       /* transverse wavenumber */
        /* x=0..1022 */
        for (x = 0; x < xmax; x += 2)
        {
            xiX = x / (xmax * 2.0 * dx) * twopi; /* transverse wavenumber */
            ga2r = k2r + xiX * xiX + xiY * xiY;
            compsqrt(ga2r, ga2i, &gar, &gai);

            bga2r = bk2r + xiX*xiX + xiY*xiY;
            compsqrt(bga2r, bga2i, &bgar, &bgai);
            compdiv(gar-bgar, gai-bgai, gar+bgar, gai+bgai, &refr, &refi);

            /* xmax=1024 -> u[0][0], u[0][1] ... u[0][1022], u[0][1023] */
            compmult(refr, refi, uin[y][x], uin[y][x+1], &urefl[y][x], &urefl[y][x+1]);
            /* xmax=1024 -> u[0][2046], u[0][2047] ... u[0][1024], u[0][1025] */
            compmult(refr, refi, uin[y][reMax-x], uin[y][imMax-x], &urefl[y][reMax-x], &urefl[y][imMax-x]);
        }
        /******* END - Handle the mid +/- frequency ******/


        /******** Handle all other frequencys ************/
        /* y = 1..511 */
        for ( y=1; y<ymax/2; y++ )
        {
            xiY = y / (ymax * dy) * twopi;

            for (x = 0; x < xmax; x += 2)
            {
                xiX = x/(xmax*2.0*dx)*twopi;
                ga2r = k2r + xiX*xiX + xiY*xiY;
                compsqrt(ga2r, ga2i, &gar, &gai);
                bga2r = bk2r + xiX*xiX + xiY*xiY;
                compsqrt(bga2r, bga2i, &bgar, &bgai);
                compdiv(gar-bgar, gai-bgai, gar+bgar, gai+bgai, &refr, &refi);
                //compmult(refr, refi, uin[y][2*x], uin[y][2*x+1], &urefl[y][2*x], &urefl[y][2*x+1]);

                /* Handle y array-index 1-511 */
                /* xmax=1024 -> u[y][0], u[y][1] ... u[y][1022], u[y][1023] */
                compmult(refr, refi, uin[y][x], uin[y][x+1], &urefl[y][x], &urefl[y][x+1]);
                /* xmax=1024 -> u[y][2046], u[y][2047] ... u[y][1024], u[y][1025] */
                compmult(refr, refi, uin[y][reMax-x], uin[y][imMax-x], &urefl[y][reMax-x], &urefl[y][imMax-x]);

                /* Handle y array-index ymax-1=1023...513. ymax/2=512 separately handled above! */
                /* xmax=1024 -> u[ymax-y][0], u[ymax-y][1] ... u[ymax-y][1022], u[ymax-y][1023] */
                compmult(refr, refi, uin[ymax-y][x], uin[ymax-y][x+1], &urefl[ymax-y][x], &urefl[ymax-y][x+1]);
                /* xmax=1024 -> u[ymax-y][2046], u[ymax-y][2047] ... u[ymax-y][1024], u[ymax-y][1025] */
                compmult(refr, refi, uin[ymax-y][reMax-x], uin[ymax-y][imMax-x], &urefl[ymax-y][reMax-x], &urefl[ymax-y][imMax-x]);
            }
        }
        /******** END - Handle all other frequencys ******/

        fourn_wrapper( urefl, xmax, ymax, true ); /* inverse spatial FFT */
        fourn_wrapper( uin, xmax, ymax, true );

        for (y=0; y<ymax; y++)
        {
            for (x=0; x<xmax; x++)
            {
                urefl[y][2*x] = refp[y][x]*urefl[y][2*x];
                urefl[y][2*x+1] = refp[y][x]*urefl[y][2*x+1];
                uin[y][2*x] = uin[y][2*x] + urefl[y][2*x];
                uin[y][2*x+1] = uin[y][2*x+1] + urefl[y][2*x+1];
            }
        }
    }
    else 
    {
        for (y=0; y<ymax; y++)
        {
            for (x=0; x<xmax; x++)
            {
                urefl[y][2*x] =   0.0; /* transmitted field, ut=T*uin = (1+R) uin */
                urefl[y][2*x+1] = 0.0;
            }
        }
    }
}


/* main program to propagate the field through the a sequence of layers */
static void propagate(sampData *sD, modelData *mD) 
{
    int x, y, z, xmax, ymax;
    double intensityTransmitted, intensityReflected, T;

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

#if ( BLOOD_DATA_TO_FILE == 1 )
    fprintfBloodData(bD, 8, 8, 8);
#endif

    xmax = sD->xAnt;
    ymax = sD->yAnt;

    /* Assign field init values.            */
    /* Input field:     Re = 1.0 Im = 0.0   */
    /* Reflected field: Re = 0.0 Im = 0.0   */
    for (y=0; y< ymax; y++)
    {
        /* E.g. sD->xAnt = 1024 s.p. x < 2047 -> 0, 2, ..., 2046 */
        for (x=0; x<(2*xmax-1); x+=2)
        {
            umig[y][x]   = initFieldValue;
            umig[y][x+1] = 0.0;
            uref[y][x]   = 0.0;
            uref[y][x+1] = 0.0;
        }
    }

    /* Calculate transmission coefficient T = 4*nrbc*nba/(nrbc+nba)^2 */
    T = (4*mD->epsilonRe*mD->backRe)/((mD->epsilonRe+mD->backRe)*(mD->epsilonRe + mD->backRe));

    intensityTransmitted = 0.0;
    /* Start propagate from z=0 to z=(SIMULATION_DEPTH-1), max{z}=(sD->zAnt-2) */
    for (z=0; z< SIMULATION_DEPTH; z++)
    {
        /* Fill two layers in xy-plane with sample points */ 
        create2DGeometry(sD, bD[z/sD->zbox], sampleLayer1, z);
        create2DGeometry(sD, bD[(z+1)/sD->zbox], sampleLayer2, z+1);

        /* Calculate initial intensity I = v*eps*E^2 = sqrt(eps)*E^2 */
        if (z == 0)
        {
            for (y=0; y<ymax; y++)
            {
                /* E.g. sD->xAnt = 1024 s.p. x=0,1,..,1023 -> max{2x+1}=2*1023 + 1 = 2047 */
                for (x=0; x<xmax; x++)
                {
                    intensityTransmitted += (umig[y][2*x] * umig[y][2*x] + umig[y][2*x+1] * umig[y][2*x+1]) *
                                             sqrt((mD->backRe + sampleLayer1[y][x] * (mD->epsilonRe - mD->backRe)) /
                                                  (mD->backRe * initFieldValue * initFieldValue));
                }
            }
        }

        /* Migrate between the two layers */
        migrate(sD, mD, sampleLayer1, umig);
        
#if ( ENABLE_BREMMER_REFLECTION == 1)
        /* Calculate interaction between two layers */
        interaction(sD, mD, sampleLayer1, sampleLayer2, umig, uref);
#endif

        /* Calculate loss in intensity due to absorption and the     */
        /* reflection. If the following condition is fullfilled      */
        /* (samplelayer1[y][x]-samplelayer2[y][x]) != 0 then there   */
        /* will be reflection between surfaces. This is valid for    */
        /* both external and internal reflection.                    */
        intensityReflected = 0.0;
        for ( y=0; y<ymax; y++)
        {
            /* E.g. sD->xAnt = 1024 s.p. x=0,1,..,1023 -> max{2x+1}=2*1023 + 1 = 2047 */
            for (x=0; x<xmax; x++)
            {
                /* Intensity loss in both external and internal reflection */
                if ((sampleLayer1[y][x] - sampleLayer2[y][x]) != 0)
                {
                    intensityReflected += (1-T) * (umig[y][2*x] * umig[y][2*x] + umig[y][2*x+1] * umig[y][2*x+1]) *
                                                   sqrt((mD->backRe + sampleLayer2[y][x] * (mD->epsilonRe - mD->backRe)) /
                                                        (mD->backRe * initFieldValue * initFieldValue));
                }
            }
        }

        /* Remove reflected intensity */
        intensityTransmitted = intensityTransmitted - intensityReflected;

#if (PRINT_MIGRATED_FIELD_TO_FILE == 1)
    if( (z % PRINT_FIELD_EVERY_NTH_LAYER) == 0)
    {
        fprintfMigratedFieldData(sD, umig, z);
    }
#endif

#if (TRANS_INTENSITY_TO_FILE == 1)
    fprintftransIntensityData(intensityTransmitted/(1.0*xmax*ymax));
#endif
        printf("z=%d\t Intensity transmitted=%.8f\n", z, intensityTransmitted/(1.0*xmax*ymax));
    }
}

/************************************************/
int main(void) 
{
    sampData sD;
    modelData mD;
    double lambda, volfrac, eps_average;

    /* Sample points for all model */	
    sD.xAnt = MODEL_TRANSVERSAL_SIZE; /* first transverse direction  */
    sD.yAnt = MODEL_TRANSVERSAL_SIZE; /* second transverse direction */
    sD.zAnt = MODEL_DEPTH_SIZE;       /* depth                       */
    sD.dx = 1;                        /* depth step                  */
    sD.dy = 1;
    sD.dz = 1;
        
    
    /* Size of box containing one bloodcell. Let this be 1/8 */
    /* of total number of samplepoints. xbox, ybox and zbox  */
    /* must be the same otherwise causing faulty geometry.   */
    /* Notice that this does not mean that xant=yant must be */
    /* equal to zant.                                        */
    sD.xbox = sD.xAnt / 8;
    sD.ybox = sD.yAnt / 8;
    sD.zbox = sD.zAnt / 8;

    /* Number of iterations for "smoothing" */
    sD.iAnt = 3;

    /* Init variables for counting occurences */
    /* in RBC vs background                   */
    nbrOfSamplespointsInRbc = 0;
    nbrOfSamplespointsInBackground = 0;
        
    /* Field data */

    /* RBC max size 7.76 um is equivalent to rbcWidthInSampPoints*/
    /* lambda = 632.8 nm = (sD.xbox/7.76)*0.6328 s.p             */

    rbcWidthInSampPoints = sD.xbox;               /* RBC width                          */
    lambda = (rbcWidthInSampPoints/7.76)*0.6328;  /* Sample points                      */
    mD.afreq = 2*pi/lambda;                       /* angular frequency c^{-1} s^{-1}    */
    mD.epsilonRe = 1.977;                         /* Re permittivity RBC                */

    /* Value chosen according to "Simulations of light scattering from a biconcave     */ 
    /* red blood cell using the finite-differencetime-domain method" by Jun Q. Lu for  */
    /* lambda = 700nm -> ni = 4.3*10^-6 -> ei = 2*nr*ni = 1.20916*10^-5                */
    mD.epsilonIm = 0.0000120916;                 /* Im permittivity */
    mD.backRe = 1.809;
    mD.backIm = 0.0;

    /* Start random number generator */
    std::srand(time(0));

    /* Model contains 2D layer with ALL samplepoints */
    propagate( &sD, &mD );

    volfrac = 1.0*nbrOfSamplespointsInRbc/(nbrOfSamplespointsInRbc + nbrOfSamplespointsInBackground);
    eps_average = (mD.epsilonRe * nbrOfSamplespointsInRbc + mD.backRe * nbrOfSamplespointsInBackground) /
        (nbrOfSamplespointsInRbc + nbrOfSamplespointsInBackground);
        
        

    /* Calculate volume fraction RBC vs background and weighted average of permitivity seen */
    printf("Volumefraction RBC vs background:%.4f percent\n", volfrac*100);
    printf("Average permitivity realpart:%.4f\n", eps_average);

    printf("Program Exit. \n");

}
/**********************************************************/
