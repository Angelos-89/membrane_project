/* the following lines go to header file later*/
double *DOS;
double dk;
int qdiag_max;
fftw_plan k2x,x2k;
void fft_setup2d(int N, double hh[], fftw_complex hq[] );
double *setup_onedspec(int N, double LL );
void free_onedspec(double S1d[]);
int wrap_around(int k,int N);
void onedspec2d(double S1d[], int N, fftw_complex hq[], double aa);
void fft();
void ifft();
/*--------header file finished ----------*/
