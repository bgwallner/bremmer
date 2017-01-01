/* Wrapper for 2D FFT calculations */

#include "fourier_ndim.h"

/* Calculate FFT */
void ddfourn(MatDoub& u1, int xmax, int ymax) {

	int j, y, x;
    int a = 0;
    
    int len = 2*xmax*ymax;
    VecDoub v1(len, a);
    
    VecInt nn(3);
    nn[0]=ymax;
    nn[1]=xmax;
    nn[2]=xmax;

	/* Copy u1[][] to v1[] */
    j = 0;
	for (y=0; y<ymax; y++) 
	{
		for (x=0; x<2*xmax; x+=2) 
		{
				v1[j] = u1[y][x] ;
				v1[j+1] = u1[y][x+1];
                //printf("Re(v1(%d)) : %f \n", j, v1[j]);
                //printf("Im(v1(%d)) : %f \n", j+1, v1[j+1]);
				j+=2;
		}
	}
	
	/* Call Inverse FFT, result returned in v2 */
	fourn(v1, nn, 1);

	/* Copy back contents from v1[] to u1[][] */
	j=0;
	for (y=0; y<ymax; y++) 
	{
		for (x=0; x<2*xmax; x+=2) 
		{
				u1[y][x] = v1[j];
				u1[y][x+1] = v1[j+1];
				j+=2;
		}
	}

	/* Deallocate temp string */

}

/* Calculate inverse FFT */
void ddfourninv(MatDoub& u1, int xmax, int ymax) {
	
	int j, y, x;
    int a = 0;
    VecDoub v1(2*xmax*ymax, a);
    VecInt nn(3);
    nn[0]=ymax;
    nn[1]=ymax;
    nn[2]=xmax;
    
	/* Copy u1[][] to v1[] */
    j = 0;
	for (y=0; y<ymax; y++) 
	{
		for (x=0; x<2*xmax; x+=2) 
		{
				v1[j] = u1[y][x];
				v1[j+1] = u1[y][x+1];
				//printf("Re(v1(%d)) : %f \n", j, v1[j]);
                //printf("Im(v1(%d)) : %f \n", j+1, v1[j+1]);
				j+=2;
		}
	}
	
	/* Call Inverse FFT, result returned in v2 */
	fourn(v1, nn, -1);

	/* Copy back contents from v1[] to u1[][] */
	j=0;
	for (y=0; y<ymax; y++) 
	{
		for (x=0; x<2*xmax; x+=2) 
		{
				u1[y][x] = v1[j];
				u1[y][x+1] = v1[j+1];
				//printf("Re(v1(%d)) : %f \n", j, v1[j]);
                //printf("Im(v1(%d)) : %f \n", j+1, v1[j+1]);
				j+=2;
		}
	}

	/* Deallocate temp string */


}