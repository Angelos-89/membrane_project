#include <valarray>
#include <fstream>
#include <fftw3.h>
#include <iostream>
#include <cmath>
#include <vector>
#include "RectMesh.hpp"

double PI = 2*acos(0.0);

/*---------------------------- FFT2d ----------------------------*/

/* This function takes as input a 2D scalar field with dimensions 
   (Nx,Ny,Nghost) and fills another RectMesh instance "Output" with 
   the values that are returned by fftw_execute. The "Output" instance 
   has dimensions (Nx+2,Ny). */

void FFT2d(RectMesh& Input_Field,RectMesh& Output)
{
  int n0 = Input_Field.getrows(); 
  int n1 = Input_Field.getcols();
  int len = n0*(n1+2);
  fftw_plan x2k;
  double* in = (double*) fftw_malloc( sizeof(double)*len );
  fftw_complex* out = (fftw_complex*) in;
  x2k = fftw_plan_dft_r2c_2d(n0,n1,in,out,FFTW_MEASURE);

  /* Fill the allocated block pointed by 
     the pointer "in" with the H data */
  
  for (int j=0; j<n0; j++){
    for(int i=0; i<n1; i++)
      in[ i + j*(n1+2) ] = Input_Field(i,j); //note the (n1+2) term 
  }

  fftw_execute(x2k);
 
  /* Fill the matrix h(q1,q2) */

  double value;
  for (int j=0; j<n0; j++){
    for(int i=0; i<n1+2; i++)
      Output(i,j) = in[i + j*(n1+2)];
  }

  /* Destroy plan & release resources */

  fftw_destroy_plan(x2k);
  fftw_free(in);
  out = nullptr;
  in  = nullptr;
}

int WrapAround(int q, int N)
{
  int k;
  if (q < N/2)
    k = q;
  else
    k = q-N;
  return k;
}

std::valarray<double> RadSpec1D(RectMesh& fft_dat,double dx)
{
  int N = fft_dat.getrows();
  double L = N*dx; //dx is the lattice spacing
  double dk = 2*PI/L; // why multiply with 2pi?
  int iq_max = sqrt(2)*(N/2); //I used N/2 beauase I take only the positive wavevectors along y axis and multiply by 2 at the end.
  std::valarray<double> ravg_spec(iq_max); //a container to store the 1D spectrum
  int k1,k2,q1,q2,re_ind,im_ind,q_sqrd,q_diag;
  double hq_re,hq_im;

  for(k1=0; k1<N/2+1; k1++) //Note that k1=0:N/2+1. I take only positive y wavelenghts
    {
      q1 = k1;
      //q1 = WrapAround(k1,N); //then this is not needed
      for(k2=0; k2<N/2+1; k2++)
	{
	  q2=k2;
	  re_ind = 2*k2;
	  im_ind = 2*k2+1;
	  hq_re = fft_dat(re_ind,k1); 
	  hq_im = fft_dat(im_ind,k1);
	  int q_sqrd = (q1*q1+q2*q2);//*(2*PI/(N*dx)); // don't remember why it should be like that. Note that this is just equal to dk, since N*dx = L.
	  q_diag = sqrt(q_sqrd);//dk; //this too.
	  if (q_diag <= iq_max)
	    ravg_spec[q_diag] += sqrt(pow(hq_re,2) + pow(hq_im,2));
	}     
    }
  return 2.0*ravg_spec/pow(N,2); // here I take into account contributions from negative wavevectors and also normalize 
}

void write2txt(std::valarray<double> data,const char* filename)
{
  std::ofstream file;
  file.open(filename);
  for(int i=0; i<data.size(); i++)
    file << data[i] << "\n";
  file.close();
}

int main()
{

  int Nx = 64;
  int Ny = 64;
  double Lx = 2*PI;
  double Ly = 2*PI;
  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double kx = 2*PI/Lx;
  double ky = 2*PI/Ly;

  std::vector<double> x,y;
  for (int i=0; i<Nx; i++)
    x.push_back(i*dx);
  for (int i=0; i<Ny; i++)
    y.push_back(i*dy);
  
  RectMesh H(Nx,Ny);
  for (int j=0; j<Ny; j++){
    for(int i=0; i<Nx; i++)
      H(i,j) =  4.0*sin(2.0*kx*x[i]) + 10.0*cos(20.0*ky*y[j]);
  }
  
  RectMesh Output(Nx+2,Ny);
  FFT2d(H,Output);
  
  std::valarray<double> rad_spec;
  rad_spec = RadSpec1D(Output,dx);

  write2txt(rad_spec,"rad_spec.txt");
  
  return 0;
}
