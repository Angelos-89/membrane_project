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
double *setup_onedspec(int N, double LL ){
  double PI =4*atan(1.);
  double *S1d;
  dk = 2*PI/LL; 
  //  int iqmax = int(sqrt(2)*N*2*PI/(N*amean))+1;
  int qdiag_max = floor(sqrt(2)*N)+1;
  S1d = (double*)malloc(qdiag_max * sizeof(double) );
  for (int iq=0;iq<qdiag_max;iq++){
    S1d[iq] = 0.;
  }
  /* put the DOS calculation here */
  return S1d;
}
/*------------------*/
void free_onedspec(double S1d[]){
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
void onedspec2d(double S1d[], int N, fftw_complex hq[], double aa ){
  int k1,k2,q1,q2,q2re,q2im;
  double hq_re,hq_im;
  /* fix the PI */
  double PI =4*atan(1.);
  for(k1=0;k1<N;k1++){
    q1 = wrap_around(k1,N);
    for(k2=0;k2<N/2+1;k2++){
      q2=k2;
      q2re = 2*k2;
      q2im = 2*k2+1;
      hq_re = hq[q2re+(N+2)*k1];
      hq_im = hq[q2im+(N+1)*k1];
      double qsqr = (2*PI/(N*aa))*(q1*q1+q2*q2); 
      int qdiag = floor(sqrt(qsqr)/dk);
      if (qdiag <= qdiag_max){
      	S1d[qdiag] = S1d[qdiag] + hq_re*hq_re+hq_im*hq_im;
      }
    }
  }
}
