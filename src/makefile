include config.angelos

objects = main.o  Vec3dLib.o RectMesh.o McmemLib.o

main : $(objects)
	$(CC) -o main $(objects) -I$(hdir_hdf5) -I$(hdir_fftw3) -L$(ldir_fftw3) $(link_fftw3) 

main.o : main.cpp Site.hpp RectMesh.hpp McmemLib.hpp
	$(CC) $(CFLAGS) main.cpp -I$(hdir_hdf5)

Vec3dLib.o : Vec3dLib.cpp Vec3dLib.hpp

RectMesh.o : RectMesh.cpp RectMesh.hpp
	$(CC) $(CFLAGS) RectMesh.cpp -I$(hdir_hdf5)

McmemLib.o : McmemLib.cpp 
	$(CC) $(CFLAGS) McmemLib.cpp -I$(hdir_hdf5)

clean :
	rm main $(objects)

