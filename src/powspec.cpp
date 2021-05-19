#include <fftw3.h>
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

  if (argc <= 1 or argc >= 3)
    {
      std::cout << "One argument is needed: Degrees of freedom per dimension. Preferably power of 2. Exiting..." << std::endl;
      exit(-1);
    }


  
  /* Generate x and y axis */
  /* ----------------------------------------------------------- */
  int N = strtol(argv[1], nullptr, 0);
  double Lx = 2*PI; //length of box
  double Ly = Lx;   
  double dx = Lx/N;
  double dy = Ly/N;
  std::vector<double> x,y;
  for (int i=0; i<N; i++)
    x.push_back(i*dx);
  y = x;
  /* ----------------------------------------------------------- */
  

  
  /* Define data for which the power spectrum will be calculated */
  /* ----------------------------------------------------------- */
  double* field = new double[N*N]{};
  std::cout << "Field is: " << std::endl;
  for (int j=0; j<N; j++)
    {
      for (int i=0; i<N; i++)
  	{
	  /* test different cases here */
  	  //field[ i + N*j ] = 1.;
	  field[ i + N*j ] = 4.0 * sin(3.0*x[i]);
	  //field[ i + N*j ] = 4.0 * cos(1.0*y[j]);
	  //field[ i + N*j ] = 5.0 * sin(-x[i] + 3.0*y[j]) + 3.0 * sin(x[i] + 3.0*y[j]);
	  //field[ i + N*j ] = 2.0 * sin( 3.0*(2*PI*x[i]/Lx) + 3.0*(2*PI*y[j]/Ly) );	    
	  std::cout << field[ i + N*j ] << "\t"; 
  	}
      std::cout << std::endl;
    }std::cout << std::endl;
  /* ----------------------------------------------------------- */



  /* Prepare in-place fft */
  /* ----------------------------------------------------------- */
  int blockSize = N*(N+2);
  fftw_plan x2k, k2x;
  double* IN = (double*) fftw_malloc( sizeof(double)*blockSize );
  fftw_complex* OUT = (fftw_complex*) IN;
  x2k = fftw_plan_dft_r2c_2d(N,N,IN,OUT,FFTW_MEASURE);
  k2x = fftw_plan_dft_c2r_2d(N,N,OUT,IN,FFTW_MEASURE);
  /* ----------------------------------------------------------- */
  


  /* Fill the allocated block pointed by the pointer "in" with the field data */
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
  double* powSpec = new double[N*(N/2+1)]{};
  double Re,Im,val,q1,q2;
  std::cout << "Power spectrum values: " << std::endl;
  for (int k2=0; k2<N; k2++){
    q2 = WrapAround(k2,N);
    for (int k1=0; k1<N/2+1; k1++){
      q1 = k1;
      Re = IN[ (N+2)*k2 + 2*k1 ];
      Im = IN[ (N+2)*k2 + 2*k1+1 ];
      val = (Re*Re + Im*Im) / pow(N,4); // normalize
      powSpec[ (N/2+1)*k2 + k1 ] = val; 
      
      // Now print the buffer
      if (val < 1e-16)
	val = 0.0;
      std::cout << "(" << q1 << "," << q2 << ")\t";
      std::cout << "(" << k1 << "," << k2 << ")\t" << val << "\t";
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
      std::cout << IN[ (N+2)*k2 + k1 ] << "\t";
    }
    std::cout << std::endl;
  }std::cout << std::endl;
  /* ----------------------------------------------------------- */



  /* Radial averaging */
  /* ----------------------------------------------------------- */
  double dk = 2*PI/Lx;
  double factor = pow(dk,2);  
  int qdiagMax = floor( sqrt( pow(dk*N/2,2) + pow(dk*N/2,2) ) ) + 1;
  int qbin;
  double qSquared;
  double* S1D = new double[qdiagMax]{}; // buffer to contain 1D spectrum
  double* DoS = new double[qdiagMax]{}; // buffer to contain the density of states
  for (int k2=0; k2<N; k2++){
    q2 = WrapAround(k2,N); 
    for (int k1=0; k1<N/2; k1++){
      q1 = WrapAround(k1,N);
      qSquared = factor * (q1*q1 + q2*q2);
      qbin = floor(sqrt(qSquared));
      S1D[qbin] += powSpec[ N/2*k2 + k1 ];
      DoS[qbin] ++;
    }
  }
  DoS[0] = 1; //the first grid point on its own bin
  DoS[qdiagMax-1] = 1; // the last grid point on its own bin
  for (int i=0; i<qdiagMax; i++)
    S1D[i] /= DoS[i];
  /* ----------------------------------------------------------- */



  /* Print results of radial averaging */
  /* ----------------------------------------------------------- */
  double printOut;
  std::cout << "Radially averaged power spectrum: " << std::endl;
  for (int i=0; i<qdiagMax; i++)
    {
      printOut = S1D[i];
      if ( S1D[i] < 1e-16 ) printOut = 0;
      std::cout << printOut << "\t";
    }std::cout << "\n" <<std::endl;

  std::cout << "Density of states: " << std::endl;
  for (int i=0; i<qdiagMax; i++)
    std::cout << DoS[i] << "\t";
  std::cout << std::endl;
  /* ----------------------------------------------------------- */

  
  
  /* Destroy plan & release resources */
  fftw_destroy_plan(x2k);
  fftw_free(IN);
  OUT = nullptr;
  IN  = nullptr;
  return 0;
}
