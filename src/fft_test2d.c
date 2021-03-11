#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
/*#include "/home/dhruba/lib/fftw3/include/fftw3.h" */
#include "fftw3.h"
#include "fft.h"
void main(){
	int N = 4;
	double *hh,*hback;
	fftw_complex *hq;
	int i,j,k, kre, kim;
//
	hh = (double*)(malloc(N*(N+2) * sizeof(double)));
//hback = (double*)(malloc(N*N * sizeof(double)));
	hq = (fftw_complex*)hh;
	fft_setup2d(N, hh, hq);
	printf(" hh in real-space \n");
	for( i=0;i<N;i++){
		for(int k=0;k<N;k++){
			hh[k+(N+2)*i] = 1.0;
			printf("%d %d %le \n",i,k,hh[k+(N+2)*i]);
		}
	}
	printf("---done\n");
	printf("in-place fft\n");
	fft();
	printf("---done\n");
	printf("hq in Fourier space (real, imag) \n");
	for(i=0;i<N;i++){
		for(k=0;k<N/2+1;k++){
			kre = 2*k;
			kim = 2*k+1;
			printf("%d %d %le %le \n",i,k,hq[kre+(N+2)*i]),hq[kim+(N+2)*i];
		}
	}
	printf("in-place inverse fft\n");
	ifft();
	printf("---done\n");
	printf(" hh back real-space (no normalization) \n");
	for(i=0;i<N;i++){
		for(k=0;k<N;k++){
			printf("%d %d %le \n",i,k,hh[k+(N+2)*i]);
		}
	}
}
/*-----------------------------------*/
