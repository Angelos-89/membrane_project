#include <fstream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
/*#include "/home/dhruba/lib/fftw3/include/fftw3.h" */
#include "RectMesh.hpp"
#include "fftw3.h"
#include "fft.h"
// int main(){
// 	int N = 4;
// 	//double *hh,*hback;
// 	fftw_complex *hq;
// 	double *hx;
// 	int i,j,k, kre, kim;
// 	RectMesh hh(N+2,N);
// 	//fftw_plan k2x,x2k; 
// 	//
// 	hx = (double*)(malloc(N*(N+2) * sizeof(double)));
// 	hq = (fftw_complex*) hx;
// 	fft_setup2d(N, hx, hq);
// 	printf(" hh in real-space \n");
// 	for( i=0;i<N;i++){
// 		for(int k=0;k<N;k++){
// 		  hh(k,i) = 1.0;
// 		  hx[k+(N+2)*i] = hh(k,i);
// 		  printf("%d %d %le \n",i,k,hh(k,i));
// 		}
// 	}
// 	printf("---done\n");
// 	printf("in-place fft\n");
// 	fft();
// 	printf("---done\n");
// 	printf("hq in Fourier space (real, imag) \n");
// 	for(i=0;i<N;i++){
// 		for(k=0;k<N/2+1;k++){
// 			kre = 2*k;
// 			kim = 2*k+1;
// 			printf("%d %d %le %le \n",i,k,
// 			      *hq[kre+(N+2)*i],*hq[kim+(N+2)*i]);
// 		}
// 	}
// 	double LL = 2*3.14159;
// 	double PI =4*atan(1.);
// 	double dk = 2*PI/LL; 
// 	int qdiag_max = floor(sqrt(2)*N)+1;
// 	printf("diag_max=%d\n",qdiag_max);
// 	double *S1d;
// 	S1d = (double*)malloc(qdiag_max * sizeof(double) );
// 	for (int iq=0;iq<qdiag_max;iq++){
// 	  S1d[iq] = 0.;
// 	}
// 	double aa = LL/N;
// 	printf(" calculate angle-summed spectra\n");
// 	onedspec2d(S1d, N, hx, aa, dk, qdiag_max);
//         printf("------done\n");
//         printf("1dspec\n");
// 	for (int iq=0; iq<qdiag_max;iq++){
// 		printf("%le \n", iq*dk);
// 	}
//         printf("------\n");
// 	printf("in-place inverse fft\n");
// 	ifft();
// 	printf("---done\n");
// 	printf(" hh back real-space (no normalization) \n");
// 	for(i=0;i<N;i++){
// 		for(k=0;k<N;k++){
// 		  printf("%d %d %le \n",i,k,hx[k+(N+2)*i] );
// 		}
// 	}
// }
/*-----------------------------------*/

double pi = 4*atan(1.0);

int main(){

  double L = 2*pi;
  int N = 16;
  int Nx = N;
  int Ny = N;
  double Lx = L;
  double Ly = L;
  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double kx = 2*pi/Lx;
  double ky = 2*pi/Ly;
  
  std::vector<double> x,y;
  for (int i=0; i<Nx; i++)
    x.push_back(i*dx);
  for (int i=0; i<Ny; i++)
    y.push_back(i*dy);
  
  RectMesh H(Nx,Ny);
  for (int j=0; j<Ny; j++){
    for(int i=0; i<Nx; i++)
      H(i,j) = 2.0*sin(2.0*kx*x[i]) + 1.0*cos(5.0*ky*y[j]);
  }


  int i,j,k, kre, kim;
  double* hx = new double[N*(N+2)]();  
  fftw_complex* hq = (fftw_complex*) hx;
  fft_setup2d(N, hx, hq);
  for( i=0;i<N;i++)
    {
      for(int k=0;k<N;k++)
	hx[k+(N+2)*i] = H(k,i);
    }
  printf("in-place fft\n");
  fft();
  printf("---done\n");
  
  double dk = 2*pi/L; 
  int qdiag_max = floor(sqrt(2)*N)+1;
  printf("diag_max=%d\n",qdiag_max);
  double *S1d = new double[qdiag_max]();
  double aa = L/(double)N;
  printf(" calculate angle-summed spectra\n");
  onedspec2d(S1d, N, hx, aa, dk, qdiag_max);
  printf("------done\n");
  
  
  std::ofstream file;
  file.open("rad_spec.txt");
  for (int i=0; i<N; i++)
    file << i*dk << "\t" << 2.0*S1d[i] << "\n";
  file.close();
  
}





