#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
/*#include "/home/dhruba/lib/fftw3/include/fftw3.h" */
#include "fftw3.h"
#define N 4 
int main(){
int n1,n2,k,i,j,kre,kim;
fftw_plan k2x,x2k;
double *omega,*omback;
fftw_complex *omk;
double *omt;

omega = (double*)(malloc(N*(N+2) * sizeof(double)));
omback = (double*)(malloc(N*N * sizeof(double)));
//omk = (fftw_complex*)(fftw_malloc(N*(N+2) * sizeof(double)));
omk = (fftw_complex*)omega;
k2x = fftw_plan_dft_c2r_2d(N,N,omk,omega,FFTW_MEASURE);
x2k = fftw_plan_dft_r2c_2d(N,N,omega,omk,FFTW_MEASURE);

for(i=0;i<N;i++){
	for(k=0;k<N;k++){
		omega[k+(N+2)*i] = 1.0;
		printf("%d %d %le \n",i,k,omega[k+(N+2)*i]);
	}
}
printf("\n");
fftw_execute(x2k);
for(i=0;i<N;i++){
	for(k=0;k<N/2+1;k++){
		kre = 2*k;
		kim = 2*k+1;
		printf("%d %d %le %le \n",i,k,omk[kre+(N+2)*i]),omk[kim+(N+2)*i];
	}
}
printf("\n");

fftw_execute(k2x);
for(i=0;i<N;i++){
	for(k=0;k<N;k++){
		omega[k+(N+2)*i] = omega[k+(N+2)*i]/((double) N*N);
	}
}

for(i=0;i<N;i++){
	for(k=0;k<N;k++){
		printf("%d %d %d %le \n",i,k,k+(N+2)*i,omega[k+(N+2)*i]);
	}
}

}
