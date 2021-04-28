include config.local
coptions := 
Include := 

ifdef FFT_IDIR
	coptions := $(coptions) -I$(FFT_IDIR)
	Include := $(Include) -I$(FFT_IDIR)
endif

ifdef HDF5_IDIR 
	coptions := $(coptions) -I$(HDF5_IDIR)
	Include := $(Include) -I$(HDF5_IDIR)
endif

ifdef FFT_LDIR 
	coptions := $(coptions) -L$(FFT_LDIR)
endif

ifdef HDF5_LDIR 
	coptions := $(coptions) -L$(HDF5_LDIR)
endif

link := -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lfftw3

#-----------------------------------------------------------------#

objects = main.o Site.o Vec3dLib.o RectMesh.o McmemLib.o fftw_lib.o

membrane.exe : $(objects)
	$(CC) -o membrane.exe $(objects) $(coptions)  $(link) 

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp $(Include) 

Site.o : Site.cpp
	$(CC) $(CFLAGS) Site.cpp 

Vec3dLib.o : Vec3dLib.cpp
	$(CC) $(CFLAGS) Vec3dLib.cpp

RectMesh.o : RectMesh.cpp
	$(CC) $(CFLAGS) RectMesh.cpp $(Include) 

McmemLib.o : McmemLib.cpp
	$(CC) $(CFLAGS) McmemLib.cpp $(Include) 

fftw_lib.o : fftw_lib.c
	$(CC) -c fftw_lib.c $(Include) 

clean :
	rm  -f membrane.exe $(objects) fftw_lib.o
