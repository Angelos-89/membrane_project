include config.local

objects = main.o Vec3dLib.o RectMesh.o McmemLib.o

membrane.exe : $(objects)
	$(CC) -o membrane.exe $(objects) -I$(HDF5.h_dir) -L$(HDF5.a_dir) -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lfftw3

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp -I$(HDF5.h_dir)

Vec3dLib.o : Vec3dLib.cpp
	$(CC) $(CFLAGS) Vec3dLib.cpp

RectMesh.o : RectMesh.cpp
	$(CC) $(CFLAGS) RectMesh.cpp -I$(HDF5.h_dir)

McmemLib.o : McmemLib.cpp
	$(CC) $(CFLAGS) McmemLib.cpp -I$(HDF5.h_dir)

clean :
	rm membrane.exe $(objects)
