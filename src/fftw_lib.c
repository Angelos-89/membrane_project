#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <fstream>
#include "fftw3.h"
#include "fft.h"

extern double PI;

int WrapAround(int k,int N)
{
  if (k <= N/2) {
    return k;}
  else{
    return k-N;}
}

void FillDoSBuffer(double DoS[],double dk,int N,int qdiagMax)
{
  double factor = pow(dk,2);  
  int q1,q2,qbin;
  double qSquared;
  for (int k2=0; k2<N; k2++)
    {
      q2 = WrapAround(k2,N); 
      for (int k1=0; k1<N/2; k1++)
	{
	  q1 = WrapAround(k1,N);
	  qSquared = factor * (q1*q1 + q2*q2);
	  qbin = floor( sqrt(qSquared)/dk );
	  DoS[qbin] ++;
	}
    }
}

void Rad1DSpec(double S1D[],double DoS[],double Hq[],double dx,double dk,int N,int qdiagMax)
{
  double LTemp = N*dx;
  double dkTemp = 2*PI / LTemp;
  double factor = pow(dkTemp,2);
  int q1,q2;
  int qBin;
  double HqRe,HqIm,qSquared;
  for(int k2=0; k2<N; k2++)
    {
      q2 = WrapAround(k2,N);
      for(int k1=0; k1<N/2+1; k1++)
  	{
  	  q1=k1;
  	  HqRe = Hq[ (N+2)*k2 + 2*k1 ] / pow(N,2);
  	  HqIm = Hq[ (N+2)*k2 + 2*k1+1 ] / pow(N,2);
  	  qSquared = factor * (q1*q1 + q2*q2);
  	  qBin = floor( sqrt(qSquared)/dk );
  	  if (qBin < qdiagMax)
  	    S1D[qBin] += pow(HqRe,2) + pow(HqIm,2);
  	}
    }
  for (int i=0; i<qdiagMax; i++)
  S1D[i] /= DoS[i];
}

void WriteSpectrum(std::string hspec_filename, double* S1D, int specSteps, int qdiagMax, double dk)
{
  std::ofstream radSpecFile;
  radSpecFile.open(hspec_filename); 
  for (int i=0; i<qdiagMax; i++)
    radSpecFile << i*dk << "\t" << S1D[i]/(double)specSteps << "\n";
  radSpecFile.close();
}




/* fftw_plan x2k; */
/* fftw_plan k2x; */
/* double dk; */
/* int qdiag_max; */
/*------------------------------------*/
/* void fft_setup2d(int N, double hh[], fftw_complex hq[]){ */
/* 	x2k = fftw_plan_dft_r2c_2d(N,N,hh,hq,FFTW_MEASURE); */
/* 	k2x = fftw_plan_dft_c2r_2d(N,N,hq,hh,FFTW_MEASURE); */
/* 	//This memory is never released with destroy plan */
/* }	 */
/* /\*---------------------*\/ */
/* void fft(){ */
/* 	fftw_execute(x2k); */
/* } */
/* void ifft(){ */
/* 	fftw_execute(k2x); */
/* } */
/* /\*---------------------*\/ */
/* void setup_onedspec(int N, double LL ){ */

/*   printf("qdiag_max=%d\n",qdiag_max); */
/* } */
/* /\*------------------*\/ */
/* void free_onedspec(double S1d[]){ */
/*   free(S1d); */
/* } */
/* wraps around */



/* // */
/* void onedspec2d(double S1d[], int N, double  hq[], double aa, double dk , int qdiag_max ) */
/* { */
/*   int k1,k2,q1,q2,q2re,q2im; */
/*   int qdiag; */
/*   double hq_re,hq_im; */
/*   // fix the PI  */
/*   for(k1=0;k1<N;k1++) */
/*     { */
/*       q1 = wrap_around(k1,N); */
/*       for(k2=0;k2<N/2+1;k2++) */
/* 	{ */
/* 	  q2=k2; */
/* 	  q2re = 2*k2; */
/* 	  q2im = 2*k2+1; */
/* 	  hq_re = hq[q2re+(N+2)*k1]/pow(N,2); */
/* 	  hq_im = hq[q2im+(N+2)*k1]/pow(N,2); */
/* 	  //	  double qsqr = (2*PI/(N*aa))*(q1*q1+q2*q2); */
/* 	  double qsqr = (4*PI*PI/(N*N*aa*aa))*(q1*q1+q2*q2); */
/* 	  //printf("h %le, %le,%le\n",qsqr, hq_re,hq_im); */
/* 	  if (qsqr == 0) */
/* 	    { */
/* 	      qdiag = 0; */
/* 	      S1d[0] += hq_re*hq_re+hq_im*hq_im; */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      qdiag = floor(sqrt(qsqr)/dk); */
/* 	      if (qdiag <= qdiag_max) */
/* 		{ */
/* 		  S1d[qdiag] += hq_re*hq_re+hq_im*hq_im; */
/* 		  //printf("qdiag, S1d %d, %le\n",qdiag, S1d[qdiag]); */
/* 		} */
/* 	    } */
	  
/* 	} */
/*     } */
/* } */

/* void onedspec2d(double S1d[], int N, double  hq[], double aa, double dk , int qdiag_max ) */
/* { */
/*   int k1,k2,q1,q2,q2re,q2im; */
/*   int qdiag; */
/*   double hq_re,hq_im; */
/*   // fix the PI  */
/*   for(k1=0;k1<N/2+1;k1++) */
/*     { */
/*       q1 = wrap_around(k1,N); */
/*       for(k2=0;k2<N/2+1;k2++) */
/* 	{ */
/* 	  q2=k2; */
/* 	  q2re = 2*k2; */
/* 	  q2im = 2*k2+1; */
/* 	  hq_re = hq[q2re+(N+2)*k1]/pow(N,2); */
/* 	  hq_im = hq[q2im+(N+2)*k1]/pow(N,2); */
/* 	  double qsqr = (2*PI/(N*aa))*(q1*q1+q2*q2); */
/* 	  //printf("h %le, %le,%le\n",qsqr, hq_re,hq_im); */
/* 	  if (qsqr == 0) */
/* 	    { */
/* 	      qdiag = 0; */
/* 	      S1d[0] += hq_re*hq_re+hq_im*hq_im; */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      qdiag = floor(sqrt(qsqr)/dk); */
/* 	      if (qdiag <= qdiag_max) */
/* 		{ */
/* 		  S1d[qdiag] += hq_re*hq_re+hq_im*hq_im; */
/* 		  //printf("qdiag, S1d %d, %le\n",qdiag, S1d[qdiag]); */
/* 		} */
/* 	    } */
	  
/* 	} */
/*     } */
/* } */

