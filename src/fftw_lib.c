#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
/*#include "/home/dhruba/lib/fftw3/include/fftw3.h" */
#include "fftw3.h"
#include "fft.h"
/*------------------------------------*/
void fft_setup2d(int N, double hh[], fftw_complex hq[] ){
	x2k = fftw_plan_dft_r2c_2d(N,N,hh,hq,FFTW_MEASURE);
	k2x = fftw_plan_dft_c2r_2d(N,N,hq,hh,FFTW_MEASURE);
}	
/*---------------------*/
void fft(){
	fftw_execute(x2k);
}
void ifft(){
	fftw_execute(k2x);
}
/*---------------------*/
void setup_onedspec(int N){
  double PI =4*atan(1.);
  //  int iqmax = int(sqrt(2)*N*2*PI/(N*amean))+1;
  int iqmax = floor(sqrt(2)*N)+1;
  S1d = (double*)malloc(iqmax * sizeof(double) );
  for (int iq=0;iq<iqmax;iq++){
    S1d[iq] = 0.;
  }
  /* put the DOS calculation here */
}
/*------------------*/
void free_onedspec(){
  free(S1d);
}
/* wraps around */
int wrap_around(int k,int N){
  if (k < N/2) {
    return k;}
  else{
    return k-N;}
  }
/*--------------------*/
void onedspec(int N, double hh[], double aa){
  int k1,k2, q2, q1,q2re,q2im;
  double hq_re,hq_im;
  /* fix the PI */
  double PI =4*atan(1.);
  for(k1=0;k1<N;k1++){
    q1 = wrap_around(k1,N);
    for(k2=0;k2<N/2+1;k2++){
      q2=k2;
      q2re = 2*k2;
      q2im = 2*k2+1;
      hq_re = hh[q2re+(N+2)*k1];
      hq_im = hh[q2im+(N+1)*k1];
      int qsqr = q1*q1+q2*q2; 
      int qdiag = floor(sqrt(qsqr)*2.*PI/(N*aa));
      S1d[qdiag] = S1d[qdiag] + hq_re*hq_re+hq_im*hq_im;
    }
  }
}
