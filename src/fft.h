#ifndef FILE_ffth_SEEN
#define FILE_ffth_SEEN
/* the following lines go to header file later*/
extern fftw_plan k2x,x2k; 
void FillDoSBuffer(double DoS[], double dk, int N, int qdiagMax);
void WriteSpectrum(std::string hspec_filename, double* S1d, int spec_steps,
		   int qdiag_max, double dk);
//void fft_setup2d(int N, double hh[], fftw_complex hq[]);
//void setup_onedspec(int N, double LL );
//void free_onedspec(double S1d[]);
int WrapAround(int k,int N);
//void onedspec2d(double S1d[], int N, double hq[], double aa, double dk, int qdiag_max);

void Rad1DSpec(double S1D[],double DoS[],double Hq[],double dx,double dk,int N,int qdiagMax);

//void fft();
//void ifft();
/*--------header file finished ----------*/
#endif /* !FILE_ffth_SEEN */
