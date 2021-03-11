/* the following lines go to header file later*/
double *S1d;
double *DOS;
fftw_plan k2x,x2k;
void fft_setup2d(int N, double hh[], fftw_complex hq[] );
void setup_onedspec(int N);
void free_onedspec();
int wrap_around(int k,int N);
void onedspec(int N, double hh[], double aa);
void fft();
void ifft();
/*--------header file finished ----------*/
