#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

double PI = 4*atan(1.0); 

int WrapAround(int k, int N)
{
  if (k <= N/2)
    return k;
  else
    return k-N;
}

int main(int argc, char* argv[])
{
  if (argc <=1) std::cout << "what is N?" << std::endl;
  int N = strtol(argv[1], nullptr, 0);
  double Lx = 2*PI;
  double Ly = Lx;
  double dx = Lx/N;
  double dy = Ly/N;
  double kx = 2*PI/Lx;
  double ky = kx;

  std::vector<double> x,y;
  for (int i=0; i<N; i++)
    x.push_back(i*dx);
  for (int i=0; i<N; i++)
    y.push_back(i*dy);

  /* Define data for which the power spectrum will be calculated */
  /* ----------------------------------------------------------- */
  double* field = new double[N*N]{};
  std::cout << "Field is: " << std::endl;
  for (int j=0; j<N; j++)
    {
      for (int i=0; i<N; i++)
  	{
  	  //field[ i + N*j ] = 1.;
	  //field[ i + N*j ] = 3.0 * sin(4.0*kx*x[i]);
	  field[ i + N*j ] = 2.0 * cos(1.0*ky*y[j]);
	  //field[ i + N*j ] = 5.0*sin(-kx*x[i] + 3.0*ky*y[j]) +
	  //                   3.0*sin(kx*x[i] + 3.0*ky*y[j]);
	  
	  std::cout << field[ i + N*j ] << "\t"; 
  	}
      std::cout << std::endl;
    }std::cout << std::endl;
  /* ----------------------------------------------------------- */

  /* Prepare in-place fft */
  /* ----------------------------------------------------------- */
  int sizeBlock = N*(N+2);
  fftw_plan x2k, k2x;
  double* IN = (double*) fftw_malloc( sizeof(double)*sizeBlock );
  fftw_complex* OUT = (fftw_complex*) IN;
  x2k = fftw_plan_dft_r2c_2d(N,N,IN,OUT,FFTW_MEASURE);
  k2x = fftw_plan_dft_c2r_2d(N,N,OUT,IN,FFTW_MEASURE);
  /* ----------------------------------------------------------- */
  
  /* Fill the allocated block pointed by the pointer "in" with the H data */
  /* ----------------------------------------------------------- */
  for (int j=0; j<N; j++){
    for(int i=0; i<N; i++)
      IN[ i + j*(N+2) ] = field[ i + N*j ]; //note the (N+2) term 
  }
  /* ----------------------------------------------------------- */

  /* Execute FFT */
  fftw_execute(x2k);

  /* Fill a buffer with the power spectrum values */
  /* ----------------------------------------------------------- */
  double* powSpec = new double[N*(N/2)]{};
  double Re,Im,val,q1,q2;
  std::cout << "Power spectrum values: " << std::endl;
  for (int k2=0; k2<N; k2++){
    q2 = WrapAround(k2,N);
    for (int k1=0; k1<N/2; k1++){
      q1 = WrapAround(k1,N);
      Re = IN[ (N+2)*k2 + 2*k1 ];
      Im = IN[ (N+2)*k2 + 2*k1+1 ];
      val = Re*Re + Im*Im;

      //-----------------------//
      if (q1==0) val = val/2.0;
      //-----------------------//

      powSpec[ k1 + (N/2)*k2 ] = val;
      
      // Now print the result
      if (val < 1e-16) val = 0.0;
      std::cout << "(" << q1 << "," << q2 << ")\t";
      std::cout << "(" << k1 << "," << k2 << ")\t" << 4*val/pow(N,4) << "\t";
    }
    std::cout << std::endl;
  }std::cout << std::endl;
  /* ----------------------------------------------------------- */

  /* Execute inverse FFT */
  fftw_execute(k2x);
  
  /* Print result of inverse FFT */
  /* ----------------------------------------------------------- */
  std::cout << "After applying inverse FFT: " << std::endl;
  for (int k2=0; k2<N; k2++){
    for (int k1=0; k1<N; k1++){
      std::cout << IN[ k1 + (N+2)*k2 ] << "\t";
    }
    std::cout << std::endl;
  }
  /* ----------------------------------------------------------- */

  /* Destroy plan & release resources */
  fftw_destroy_plan(x2k);
  fftw_free(IN);
  OUT = nullptr;
  IN  = nullptr;
  return 0;
}
