#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
/*#include "/home/dhruba/lib/fftw3/include/fftw3.h" */
#include "fftw3.h"
#include "fft.h"
fftw_plan x2k;
fftw_plan k2x;
double dk;
int qdiag_max;
extern double PI;
/*------------------------------------*/
void fft_setup2d(int N, double hh[], fftw_complex hq[]){
	x2k = fftw_plan_dft_r2c_2d(N,N,hh,hq,FFTW_MEASURE);
	k2x = fftw_plan_dft_c2r_2d(N,N,hq,hh,FFTW_MEASURE);
	//This memory is never released with destroy plan
}	
/*---------------------*/
void fft(){
	fftw_execute(x2k);
}
void ifft(){
	fftw_execute(k2x);
}
/*---------------------*/
void setup_onedspec(int N, double LL ){

  printf("qdiag_max=%d\n",qdiag_max);
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
void onedspec2d(double S1d[], int N, double  hq[], double aa, double dk , int qdiag_max )
{
  int k1,k2,q1,q2,q2re,q2im;
  int qdiag;
  double hq_re,hq_im;
  // fix the PI 
  for(k1=0;k1<N;k1++)
    {
      q1 = wrap_around(k1,N);
      for(k2=0;k2<N/2+1;k2++)
	{
	  q2=k2;
	  q2re = 2*k2;
	  q2im = 2*k2+1;
	  hq_re = hq[q2re+(N+2)*k1]/pow(N,2);
	  hq_im = hq[q2im+(N+2)*k1]/pow(N,2);
	  double qsqr = (2*PI/(N*aa))*(q1*q1+q2*q2);
	  //printf("h %le, %le,%le\n",qsqr, hq_re,hq_im);
	  if (qsqr == 0)
	    {
	      qdiag = 0;
	      S1d[0] += hq_re*hq_re+hq_im*hq_im;
	    }
	  else
	    {
	      qdiag = floor(sqrt(qsqr)/dk);
	      if (qdiag <= qdiag_max)
		{
		  S1d[qdiag] += hq_re*hq_re+hq_im*hq_im;
		  //printf("qdiag, S1d %d, %le\n",qdiag, S1d[qdiag]);
		}
	    }
	  
	}
    }
}

