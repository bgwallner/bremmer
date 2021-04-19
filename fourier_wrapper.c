/* Wrapper for 2D FFT calculations. Field is represented as */
/* MatDoub[][] while fourn() input expect VecDoub[].        */

#include "fourier_ndim.h"
#include "stdlib.h"

/* Calculate FFT */
void fourn_wrapper( MatDoub& u1, int xmax, int ymax, bool inverse ) 
{

	int i, y, x;

	/* Create a Vector with double and all zeros */
	VecDoub v1(2*xmax*ymax, 0.0);
    
	/* 2D input to NR fourn() */
    VecInt nn(2);
    nn[0]=ymax;
    nn[1]=2*xmax;

	/* Copy u1[ymax][2*xmax] to v1[2*xmax*ymax] */
    i = 0;
	for (y=0; y<ymax; y++) 
	{
		for (x=0; x<2*xmax; x+=2) 
		{
				v1[i] = u1[y][x] ;
				v1[i+1] = u1[y][x+1];
				i+=2;
		}
	}
	
    if ( false == inverse )
    {
        /* Call FFT */
        fourn(v1, nn, 1);
	}
    else
    {
        /* Inverse FFT */
        fourn(v1, nn, -1);
    }

	/* Copy back contents from v1[2*xmax*ymax] to u1[ymax][2*xmax] */
	i=0;
	for (y=0; y<ymax; y++) 
	{
		for (x=0; x<2*xmax; x+=2) 
		{
				u1[y][x] = v1[i];
				u1[y][x+1] = v1[i+1];
				i+=2;
		}
	}
}
