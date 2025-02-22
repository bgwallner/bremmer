/* Bo-Göran Wallner 2021 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <random>
#include <complex>
#include "nr3.h"
#include "fourier_ndim.h"

/* Some constants*/
#define twopi                  6.283185307177959
#define pi                     3.14159265358
#define max32bit               4294967295
#define initFieldValue         1.0
#define twoPiInDegrees         360

/* Development environment */
#define GNU_LINUX              0

/******************* Optical properties *********************/
/* Reference values according to "Simulations of light scattering from a biconcave */ 
/* red blood cell using the finite-differencetime-domain method" by Jun Q. Lu for  */
/* lambda = 700nm -> ni = 4.3*10^-6 -> ei = 2*nr*ni = 1.20916*10^-5                */
/* lambda = 1000nm -> ni = 1.68*10^-5 -> ei = 6.64272*10^-5                        */

/* Real-part values are chosen from "Numerical simulations of light scattering by  */
/* red blood cells" by Karlsson, Anders etc.                                       */

/* Absorption values chosen according to article:                                  */
/* "A literature review and novel theoretical approach on the optical properties   */
/*  of whole blood" by Nienke Bosschaart etc. This article also provide scattering */
/* coefficient to be able to compare with the Bremmer series and to use consistent */
/* values.                                                                         */
/* lambda = 635nm  -> ua = 0.24 mm^-1 -> ei = 3.3984678*10^-5                      */
#define LIGHT_WAVE_LENGTH      0.6328  /* um */
#define RBC_WIDTH              7.76    /* um */
#define RBC_PERMITIVITY_RE     1.977
#define RBC_PERMITIVITY_IM     0.000033984678
#define BA_PERMITIVITY_RE      1.809
#define BA_PERMITIVITY_IM      0.0

/* RBC geometry parametric values */
#define _R0  3.88
#define _C0  0.81
#define _C2  7.83
#define _C4 -4.39

/* Use flag for simulating a one-dimensional slab. Possible */
/* to compare with Fresnels equations for testing method.   */
#define ONE_DIMENSIONAL_SLAB 0
#define RECTANGLE_START      256
#define RECTANGLE_END        768
#define SLAB_START           512
#define SLAB_END             768

/* Use flag to have spherical symmetry instead of RBC shape */
#define GEOMETRY_SPHERICAL 0
/* Radius of sphere giving eqv RBC volume 94 um^2 */
#define SPHERE_RADIUS_EQ   46.52 /* radius=2.82 um and MODEL_DIMENSION=1024 */

/* These defines can be changed to modify the number of sampling-points     */
/* of the whole model and the number of RBCs. Care must be taken so that    */
/* the wavelength is sampled 10 times or more. E.g. MODEL_DIMENSION =       */
/* 1024 and NBR_RBC_IN_ONE_ROW = 8 will give wavelength ~10 samplingpoints. */
/* If choosing MODEL_DIMENSION=1024 and having NBR_RBC_IN_ONE_ROW = 16      */
/* would cause wavelength to be ~5 samplingpoints which is to small.        */

/* Model size chosen as 2^N*1024, N=0,1,2...        */
/* NOTICE! Log-files have a huge disk-demand and    */
/* if disk runs out of space the program segfaults. */
/* Using 1024 roughly give 10GB of data.            */
/* Verified to work up to N=3 -> 8192x8192 matrix   */
/* while N=4 causes segmentation fault. May be      */
/* due to fundamental datatype overflowing.         */
#define MODEL_DIMENSION        1024
/* Number of RBCs in one row, 2^N, N=0,1,2...       */
#define NBR_RBC_IN_ONE_ROW     8

/* Number of iterations in migrate() for minimizing change of    */
/* wavenumber due to difference in realpart of permitivity.      */
/* Each calculation step is calculated eps0+d, eps0+2d,...,eps1  */
/* get a smooth transition from eps0->eps1 in steps of some d.   */
#define NBR_SMOOTING_ITERATIONS 5

/* Choose *one* only angle setup for different simulations */
#define ANGLES_RANDOM          0
#define ANGLE_THETA_ZERO       1
#define ANGLE_THETA_PI_HALF    0
#define ANGLE_THETA_PI_FOURTH  0

/* When using fix theta the angles psi and fi can */
/* either be random or fix 0.                     */
#define ANGLE_RANDOM_PSI_FI    1

/* Enable backward field */
#define ENABLE_BREMMER_REFLECTION 1

/* Enable random xy-plane translation of RBCs */
#define ENABLE_RANDOM_TRANSLATION 1

/* Enable custom RBC width. Otherwise is automatically */
/*  set to max size sD.xAnt / NBR_RBC_IN_ONE_ROW       */
#define ENABLE_RBC_WIDTH_CUSTOM   0
/* Custom width in sampling points */
#define RBC_CUSTOM_WIDTH          256

/* To be able som simulate depths larger than MODEL_DIMENSION */
/* can run consequtive simulation using field output from     */
/* previous simulation as initial value for the field. For    */
/* each new simulation a new random 3D geometrical model is   */
/* created. Arbitrary number of simulations can be created.   */
/* Care should be taken for disk-space if logging fields.     */
/* Values shall be 1,2,...N.                                  */
#define MAX_MODEL_NUMBER 1

/* Number of layers to simulate in propagation direction */
#define SIMULATION_DEPTH       MODEL_DIMENSION-1

/**** Flags for debug output in textfiles ****/

/* Print multipe geometrical layer of RBC to file. If vertical  */
/* field shall be plotted against geometrical model the value   */
/* shall be set equal to SIMULATION_DEPTH to be able to plot    */
/* in e.g. MATLAB. E.g. MODEL_DIMENSION=1024 will give          */
/* output samplelayer1.txt -> samplelayer1022.txt.              */
#define PRINT_GEOMETRY_TO_FILE       1
/* Print geometry between START and STOP */
#define START_PRINT_GEOMETRY_TO_FILE 0
#define STOP_PRINT_GEOMETRY_TO_FILE  SIMULATION_DEPTH

/* Print all centers, angles etc for each RBC */
#define BLOOD_DATA_TO_FILE                  0

/* Print transmitted intensity I=I(z) */
#define TRANS_INTENSITY_TO_FILE             1

/* Print electric field to file, real-part and imaginary part */
#define PRINT_MIGRATED_FIELD_TO_FILE        1
#define PRINT_FIELD_START                   0
#define PRINT_FIELD_STOP                    SIMULATION_DEPTH
/* Using this flag will print sqrt(eps)*E^2 instead of real-part */
#define PRINT_MIGRATED_INTENSITY            1

/* Set in which model to print fields in   */
/* Default set to "last" simulation model. */
#define MODEL_NBR_TO_PRINT_FIELD_IN         MAX_MODEL_NUMBER-1

/* File handles */
static FILE *fpSampleLayer;
static FILE *fpBloodData;
static FILE* fpTransPower;
static FILE* fpInteractionData;
static FILE* fpMigrateData;
static FILE* fpMigratedField;
static FILE* fpImaginaryField;

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
static void fprintfMigratedFieldData(sampData* sD, modelData* mD, MatDoub& field, MatInt& sampleLayer, int z)
{
    int x, y;
    double sqrtEps2;
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
#if ( PRINT_MIGRATED_INTENSITY == 1 )
            /* Print transmitted intensity distribution */
            sqrtEps2 = sqrt((mD->backRe + sampleLayer[y][x] * (mD->epsilonRe - mD->backRe)) / mD->backRe);
            fprintf(fpMigratedField, "%.4f\t", sqrtEps2 *(field[y][2*x]*field[y][2*x] + field[y][2*x+1]*field[y][2*x+1]) /(initFieldValue*initFieldValue));
#else
            /* Print real-part of the field */
            fprintf(fpMigratedField, "%.4f\t", field[y][2*x]/initFieldValue);
#endif
        }
        fprintf(fpMigratedField, "\n");
    }
    fclose(fpMigratedField);
}

/* Function prints one migrated field data to file */
static void fprintfImaginaryFieldData(sampData* sD, MatDoub& field, int z)
{
    int x, y;
    char buf[40];
    snprintf(buf, sizeof(buf), "data/fieldData/imaginaryfield%d.txt", z);

#if ( GNU_LINUX == 1 )
    fpImaginaryField = fopen(buf, "w");
#else
    fopen_s(&fpImaginaryField, buf, "w");
#endif

    for (y = 0; y < sD->yAnt; y++)
    {
        for (x = 0; x < sD->xAnt; x++)
        {
            /* Print real-part of the field */
            fprintf(fpImaginaryField, "%.4f\t", field[y][2*x+1]/initFieldValue);
        }
        fprintf(fpImaginaryField, "\n");
    }
    fclose(fpImaginaryField);
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
    for (i=0; i<dim1; i++)
    {
        a[i] = (bloodData ***) malloc(dim2 * sizeof(bloodData **));
        for (j=0; j<dim2; j++)
        {
            a[i][j] = (bloodData **) malloc(dim3 * sizeof(bloodData *));
            for (k=0; k<dim3; k++)
            {
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

    dx = 0;
    dy = 0;
    for (zb = 0; zb < nbrOfZbox; zb++)
    {
        /* Displace centers of whole xy-plane layer         */
        /* in y-direction to create randomness in geometry. */
#if (ENABLE_RANDOM_TRANSLATION == 1)
        dy = (int)random1((double)sD->ybox);
#endif
        for (yb = 0; yb < nbrOfYbox; yb++)
        {
            /* For every fix yb translate in x-direction */
#if (ENABLE_RANDOM_TRANSLATION == 1)
            dx = (int)random1((double)sD->xbox);
#endif
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

                /* Init values */
                bD[zb][yb][xb]->theta = 0.0;
                bD[zb][yb][xb]->fi    = 0.0;
                bD[zb][yb][xb]->psi   = 0.0;
#if    (ANGLES_RANDOM == 1)
                bD[zb][yb][xb]->theta = random1(pi);
                bD[zb][yb][xb]->fi = random1(2*pi);
                bD[zb][yb][xb]->psi = random1(2*pi);
#endif

#if    (ANGLE_THETA_ZERO == 1)
                bD[zb][yb][xb]->theta = 0.0;
#if    (ANGLE_RANDOM_PSI_FI == 1)
                bD[zb][yb][xb]->fi = random1(2*pi);
                bD[zb][yb][xb]->psi = random1(2*pi);
#endif
#endif

#if    (ANGLE_THETA_PI_HALF == 1)
                bD[zb][yb][xb]->theta = pi/2;
#if    (ANGLE_RANDOM_PSI_FI == 1)
                bD[zb][yb][xb]->fi = random1(2*pi);
                bD[zb][yb][xb]->psi = random1(2*pi);
#endif
#endif

#if    (ANGLE_THETA_PI_FOURTH == 1)
                bD[zb][yb][xb]->theta = pi/4;
#if    (ANGLE_RANDOM_PSI_FI == 1)
                bD[zb][yb][xb]->fi = random1(2*pi);
                bD[zb][yb][xb]->psi = random1(2*pi);
#endif
#endif
            }
        }
    }
}

/* 
* Fill 2D geometry layer fox fixed z 
*/

static void create2DGeometry(sampData *sD, bloodData ***bD, MatInt& geometry2D, int z, int modelNumber)
{
    /* Coordinates in global coordinate system */
    int x, y, zb, xb, dx, dy, yb, cellsize, xmax, ymax;
    /* Coordinates in local cell with disk in origo, e.g. xl = f(x,y,z) */
    int xl, yl, zl;
    /* Coordinates in Euler rotated system, e.g. xe = f(xl,yl,zl)       */
    int xe, ye, ze, Re;
    /* Blood-cell constants */
    double C0, C2, C4, R0, D_xDiff;
    /* Trigonometric values */
    double cost, sint, sinf, cosf, sinp, cosp;

    xmax = sD->xAnt;
    ymax = sD->yAnt;
    cellsize = sD->xbox;

    /* Iterate over the whole xy-plane */
    for (y=0; y<ymax; y++)
    {
        for (x=0; x<xmax; x++)
        {
            /* Calculate which cell */
            zb=z/cellsize;
            xb=x/cellsize;
            yb=y/cellsize;

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

            R0 = _R0*rbcWidthInSampPoints/RBC_WIDTH;
            C0 = _C0*rbcWidthInSampPoints/RBC_WIDTH;
            C2 = _C2*rbcWidthInSampPoints/RBC_WIDTH;
            C4 = _C4*rbcWidthInSampPoints/RBC_WIDTH;
            /* Calculate (xe,ye) corresponding radius */
            Re = sqrt(pow(xe,2) + pow(ye,2));

            /* The whole xy-layer shall be translated */
            dx = bD[yb][xb]->dx;
            dy = bD[yb][xb]->dy;

#if (ONE_DIMENSIONAL_SLAB == 1)
                /* One dimensional slab */
                if ( z <= (SLAB_END) && (z >= SLAB_START) )
                {
                    if (x >= RECTANGLE_START && x <= RECTANGLE_END && 
                        y >= RECTANGLE_START && y <= RECTANGLE_END)
                    {
                        geometry2D[x][y] = 1;
                        nbrOfSamplespointsInRbc++;

                    }
                    else
                    {
                        geometry2D[x][y] = 0;
                        nbrOfSamplespointsInBackground++;
                    }
                    //geometry2D[x][y] = 1;
                    //nbrOfSamplespointsInRbc++;
                }
                else
                {
                    geometry2D[x][y] = 0;
                    nbrOfSamplespointsInBackground++;
                }
#else
            /* Check if (xe,ye) is within disk since expression only valid within */
            /* it and that exclude e.g. corners of the cell.                      */
            if(Re < rbcWidthInSampPoints/2)
            {
#if ( GEOMETRY_SPHERICAL == 1 )
                /* Spherical symmetry */
                D_xDiff = pow(xe,2) + pow(ye,2) + pow(ze,2) - pow(SPHERE_RADIUS_EQ,2);
#else
                /* RBC - Re is valid to use in the expression */
                D_xDiff = 4*pow(ze,2) - (1-pow(Re/R0,2)) * pow(C0 + C2*pow(Re/R0,2) + C4*pow(Re/R0,4), 2);
#endif
                /* Check if ze is outside or within the biconcave disk */
                if ( D_xDiff > 0.0)
                {
                    /* ze is outside the disk */
                    geometry2D[(x+dx)%xmax][(y+dy)%ymax] = 0;
                    nbrOfSamplespointsInBackground++;
                }
                else
                {
                /* ze is within the biconcave disk */
                    geometry2D[(x+dx)%xmax][(y+dy)%ymax] = 1;
                    nbrOfSamplespointsInRbc++;
                }
            }
            else
            {
               /* Re is outside the disk */
                geometry2D[(x+dx)%xmax][(y+dy)%ymax] = 0;
                nbrOfSamplespointsInBackground++;
            }
#endif
        }
    }

#if ( PRINT_GEOMETRY_TO_FILE == 1)
    if ( modelNumber == MODEL_NBR_TO_PRINT_FIELD_IN )
    {
        /* Print layers with to files */
        if ( (z >= START_PRINT_GEOMETRY_TO_FILE) && 
             (z <= STOP_PRINT_GEOMETRY_TO_FILE) )
        {
            fprintfsampData(sD, geometry2D, z);
        }
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
        epsi = -s* mD->epsilonIm;                                             /* Im permittivity*/
        compmult(eta, w, eta, w, &kr, &ki);                                   /* wave number, c^{-1} s */
        compmult(kr, ki, epsr, epsi, &k2r, &ga2i);                            /* sqr wave number, (kr+i*ki)(eps + i*epsi), .^{2} */

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

#if ( ENABLE_BREMMER_REFLECTION == 1)
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
#endif /* ENABLE_BREMMER_REFLECTION */
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
    int x, y, z, xmax, ymax, modelNbr, i, samples;
    double intensityTransmitted, intensityReflected, intensityTimeTrans, intensityTimeRefl;
    double EphasorReT, EphasorImT, EphasorReR, EphasorImR, wt, sqrtEps2, sqrtEps1, sin2wt, cos2wt;

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

    /* Numbers of models to simulate */
    for (modelNbr=0; modelNbr<MAX_MODEL_NUMBER; modelNbr++ )
    {
        /* Initialize the 3D grid, blood cell centers only */
        initBloodData3DArray( sD, bD );

        #if ( BLOOD_DATA_TO_FILE == 1 )
            fprintfBloodData(bD, NBR_RBC_IN_ONE_ROW, NBR_RBC_IN_ONE_ROW, NBR_RBC_IN_ONE_ROW);
        #endif /* BLOOD_DATA_TO_FILE */

        /* Start propagate from z=0 to z=(SIMULATION_DEPTH-1) */
        for (z=0; z<SIMULATION_DEPTH; z++)
        {
            /* Fill two layers in xy-plane with sample points */ 
            create2DGeometry(sD, bD[z/sD->zbox], sampleLayer1, z, modelNbr);
            create2DGeometry(sD, bD[(z+1)/sD->zbox], sampleLayer2, z+1, modelNbr);

            /* Migrate between the two layers */
            migrate(sD, mD, sampleLayer1, umig);

            /* Calculate interaction between two layers */
            interaction(sD, mD, sampleLayer1, sampleLayer2, umig, uref);

            intensityReflected   = 0.0;
            intensityTransmitted = 0.0;

            /* Calculate intensity (time-averaged Poynting vector) from phasors.  */
            for ( y=0; y<ymax; y++)
            {
                /* E.g. sD->xAnt = 1024 s.p. x=0,1,..,1023 -> max{2x+1}=2*1023 + 1 = 2047 */
                for (x=0; x<xmax; x++)
                {
                    sqrtEps2 = sqrt((mD->backRe + sampleLayer2[y][x] * (mD->epsilonRe - mD->backRe)) / mD->backRe);
                    sqrtEps1 = sqrt((mD->backRe + sampleLayer1[y][x] * (mD->epsilonRe - mD->backRe)) / mD->backRe);

                    /* Time-dependent physical field E(t)=E'cos(wt)+E''sin(wt) */
                    EphasorReT = umig[y][2*x];   /* Transmitted E'  */
                    EphasorImT = umig[y][2*x+1]; /* Transmitted E'' */
                    EphasorReR = uref[y][2*x];   /* Reflected E'    */
                    EphasorImR = uref[y][2*x+1]; /* Reflected E''   */

                    /* Intensity = n*c*E(t)^2 ~ sqrt(eps)*(E'cos(wt))^2 + (E''sin(wt))^2 + E'E''sin(2wt) */
                    intensityTimeTrans = 0.0;
                    intensityTimeRefl  = 0.0;
 
                    /* <Intensity> = nc x EE* ~ sqrt(eps)*(E'^2 + E''^2) */
                    intensityTimeTrans = (EphasorReT * EphasorReT + EphasorImT * EphasorImT) * sqrtEps2;
                    intensityTimeRefl = (EphasorReR * EphasorReR + EphasorImR * EphasorImR) * sqrtEps1;

                    /* Add to intensities for whole layer */
                    intensityTransmitted += intensityTimeTrans;
                    intensityReflected += intensityTimeRefl;
                }
            }

    #if (PRINT_MIGRATED_FIELD_TO_FILE == 1)
        if ( modelNbr == MODEL_NBR_TO_PRINT_FIELD_IN )
        {
            if( (z >= PRINT_FIELD_START) && (z <= PRINT_FIELD_STOP) )
            {
                fprintfMigratedFieldData(sD, mD, umig, sampleLayer2, z);
                fprintfImaginaryFieldData(sD, umig, z);
            }
        }
    #endif /* PRINT_MIGRATED_FIELD_TO_FILE */

    #if (TRANS_INTENSITY_TO_FILE == 1)
        fprintftransIntensityData(intensityTransmitted/(1.0*xmax*ymax));
    #endif /* TRANS_INTENSITY_TO_FILE */
            printf("z=%d\t Transmitted powerflux:      Pt=%.8f\n", (z+1+modelNbr*MODEL_DIMENSION), intensityTransmitted/(1.0*xmax*ymax));
            printf("\t Reflected powerflux:        Pr=%.8f\n", intensityReflected/(1.0*xmax*ymax));
        }
    }
}

/************************************************/
int main(void) 
{
    sampData sD;
    modelData mD;
    time_t startTime, stopTime;
    double lambda, volfrac, eps_average;

    /* Benchmarking */
    startTime = time(NULL);
    stopTime  = time(NULL);
    //ctime(&startTime);

    /* Sample points for all model. Always equal size.        */	
    sD.xAnt = MODEL_DIMENSION; /* first transverse direction  */
    sD.yAnt = MODEL_DIMENSION; /* second transverse direction */
    sD.zAnt = MODEL_DIMENSION; /* depth                       */
    sD.dx = 1;                 /* depth step                  */
    sD.dy = 1;
    sD.dz = 1;
        
    
    /* Size of box containing one bloodcell. Let this be 1/NBR_RBC_IN_ONE_ROW  */
    /* of total number of samplepoints. Always equal size.                     */
    sD.xbox = sD.xAnt / NBR_RBC_IN_ONE_ROW;
    sD.ybox = sD.yAnt / NBR_RBC_IN_ONE_ROW;
    sD.zbox = sD.zAnt / NBR_RBC_IN_ONE_ROW;

    /* Number of iterations for "smoothing" */
    sD.iAnt = NBR_SMOOTING_ITERATIONS;

    /* Init variables for counting occurences */
    /* in RBC vs background                   */
    nbrOfSamplespointsInRbc = 0;
    nbrOfSamplespointsInBackground = 0;

#if ( ENABLE_RBC_WIDTH_CUSTOM == 1)
    rbcWidthInSampPoints = RBC_CUSTOM_WIDTH;                      /* RBC custom width */
#else
    rbcWidthInSampPoints = sD.xbox;                               /* RBC max width    */
#endif /* ENABLE_RBC_WIDTH_CUSTOM */

    /******************* Field data ******************/
    lambda = (rbcWidthInSampPoints/RBC_WIDTH)*LIGHT_WAVE_LENGTH;  /* Sample points                      */
    mD.afreq = 2*pi/lambda;                                       /* angular frequency c^{-1} s^{-1}    */
    mD.epsilonRe = RBC_PERMITIVITY_RE;                            /* Re permittivity RBC                */
    mD.epsilonIm = RBC_PERMITIVITY_IM;                            /* Im permittivity                    */
    mD.backRe = BA_PERMITIVITY_RE;
    mD.backIm = BA_PERMITIVITY_IM;

    /* Start random number generator */
    std::srand(time(0));

    /* Start propagation of wave */
    propagate( &sD, &mD );

    volfrac = 1.0*nbrOfSamplespointsInRbc/(nbrOfSamplespointsInRbc + nbrOfSamplespointsInBackground);
    eps_average = (mD.epsilonRe * nbrOfSamplespointsInRbc + mD.backRe * nbrOfSamplespointsInBackground) /
        (nbrOfSamplespointsInRbc + nbrOfSamplespointsInBackground);

    stopTime = time(NULL);
    /* Calculate volume fraction RBC vs background and weighted average of permitivity seen */
    printf("Volumefraction RBC vs background:%.4f percent\n", volfrac*100);
    printf("Average permitivity realpart:%.4f\n", eps_average);
    printf("Wavelength in samplingpoints:%.8f\n", lambda);
    printf("Angular frequency:%.8f\n", mD.afreq);
    printf("######### BENCHMARKING #########\n");
    //printf("Simulation started: %s", ctime(&startTime));
    //printf("Simulation ended: %s", ctime(&stopTime));
    printf("Program Exit. \n");

}
/**********************************************************/
