include config.local
#FFT_idir = /home/dhruba/Library/include/
#FFT_ldir = /home/dhruba/Library/lib/

objects = main.o Site.o Vec3dLib.o RectMesh.o McmemLib.o fftw_lib.o
ffto = test_2dfft.o fftw_lib.o RectMesh.o  

# r_spec.exe : $(objects)
# 	$(CC) -o r_spec.exe $(objects) -I$(HDF5.h_dir) -L$(HDF5.a_dir) -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lfftw3

#ftest : $(ffto)
#	gcc -o ftest.exe $(ffto) -l fftw3  -I$(FFT_idir) -L$(FFT_ldir) -lm 

membrane.exe : $(objects)
	$(CC) -o membrane.exe $(objects) -I$(HDF5.h_dir) -L$(HDF5.a_dir) -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lfftw3

ftest : $(ffto)
	$(CC) -o ftest.exe $(ffto) -I$(HDF5.h_dir) -L$(HDF5.a_dir) -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lfftw3

r_spec.o : r_spec.cpp
	$(CC) $(CFLAGS) r_spec.cpp -I$(HDF5.h_dir)

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp -I$(HDF5.h_dir)

Site.o : Site.cpp
	$(CC) $(CFLAGS) Site.cpp 

Vec3dLib.o : Vec3dLib.cpp
	$(CC) $(CFLAGS) Vec3dLib.cpp

RectMesh.o : RectMesh.cpp
	$(CC) $(CFLAGS) RectMesh.cpp -I$(HDF5.h_dir)

McmemLib.o : McmemLib.cpp
	$(CC) $(CFLAGS) McmemLib.cpp -I$(HDF5.h_dir)

fftw_lib.o : fftw_lib.c
	$(CC) -c fftw_lib.c -I$(FFT_idir)

fft_test2d.o : fft_test2d.c
	gcc -c fft_test2d.c -I$(FFT_idir)

test_2dfft.o : test_2dfft.cpp
	$(CC) -c test_2dfft.cpp -I$(FFT_idir)

clean :
	rm  -f membrane.exe ftest.exe $(objects) $(ffto)
