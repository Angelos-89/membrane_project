#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "RectMesh.hpp"

/*--------------------- Constructor -----------------------------*/

void RectMesh::CheckConstr(int& cols, int& rows, int& ghost)
{
  if (cols < 2 || rows < 2 || ghost < 0)
    {
      std::cout << "The smallest mesh must have dimensions (2x2) "
	"and the minimum number of ghost points is zero. Exiting."
		<< std::endl; 
      exit(EXIT_FAILURE);
    }
}

RectMesh::RectMesh(int cols, int rows, int ghost)
  : cols_x(cols), rows_y(rows), nghost(ghost), mesh(nullptr)
{
  CheckConstr(cols,rows,ghost);
  int len = (cols_x+2*nghost)*(rows_y + 2*nghost);
  mesh = (double*) fftw_malloc(len*sizeof(double));
  if (mesh != nullptr){
    std::fill(mesh, mesh+len, 0.0);
    //std::cout << "RectMesh constructor called." << std::endl;
  }
  else
    {
      std::cout << "RectMesh allocation failed. Exiting. " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

/*------------------- Destructor -----------------------*/

RectMesh::~RectMesh()
{
  cols_x = rows_y = nghost = 0;
  fftw_free(mesh);
  mesh = nullptr;
  //std::cout << "RectMesh destructor called." << std::endl;
}

/*----------------- Copy constructor ---------------------------*/

RectMesh::RectMesh(const RectMesh& obj)
  : cols_x(obj.cols_x), rows_y(obj.rows_y),
    nghost(obj.nghost), mesh(nullptr) 
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  mesh = (double*) fftw_malloc(len*sizeof(double));
  std::copy(obj.mesh, obj.mesh+len, mesh);
  //std::cout << "RectMesh copy constructor called. " << std:: endl;
}

/*------------- Move Constructor ---------------*/

RectMesh::RectMesh(RectMesh&& tmp) noexcept
  : cols_x(tmp.cols_x), rows_y(tmp.rows_y),
    nghost(tmp.nghost), mesh(tmp.mesh)
{
  tmp.cols_x = tmp.rows_y = tmp.nghost = 0;
  tmp.mesh = nullptr;
  //std::cout << "RectMesh move constructor called." << std::endl;
}

/*--------- Copy assignment operator -----------*/

RectMesh& RectMesh::operator=(const RectMesh& rhs)
{ 
  if (this != &rhs)
    {
      int len = (rhs.cols_x+2*rhs.nghost)*(rhs.rows_y+2*rhs.nghost);
      if (mesh!= nullptr) fftw_free(mesh);
      cols_x = rhs.cols_x;
      rows_y = rhs.rows_y;
      nghost = rhs.nghost;
      mesh = (double*) fftw_malloc(len*sizeof(double));
      std::copy(rhs.mesh, rhs.mesh+len, mesh);
    }
  //std::cout << "RectMesh copy assignment operator called." << std::endl;
  return *this;
}

/*---------------- Move assignment operator -------------*/

RectMesh& RectMesh::operator=(RectMesh&& tmp)
{
  if (this != &tmp)
    {
      cols_x = tmp.cols_x;
      rows_y = tmp.rows_y;
      nghost = tmp.nghost;
      mesh = tmp.mesh;
      tmp.cols_x = tmp.rows_y = tmp.nghost = 0;
      tmp.mesh = nullptr;
    }
  //std::cout << "RectMesh move assignment operator called." << std::endl;
  return *this;
}

/*------------- Operator () overloading ---------*/

double& RectMesh::operator()(int i, int j) const
{
  int jump = cols_x + 2*nghost;
  return mesh[ (i+nghost) + (j+nghost)*jump ];
}

double& RectMesh::operator()(int i, int j)
{
  int jump = cols_x + 2*nghost;
  return mesh[ (i+nghost) + (j+nghost)*jump ];
}

/*-------------------- Algebra ------------------------*/

RectMesh RectMesh::operator+(const RectMesh& rhs) const
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  RectMesh temp(rows_y+2*nghost, cols_x+2*nghost);
  for (int i=0; i<len; i++)
    temp.mesh[i] = mesh[i] + rhs.mesh[i];
  return temp;
}


RectMesh& RectMesh::operator+=(const RectMesh& rhs)
{
  //check if dimensions agree?
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  for(int i=0; i<len; i++)
    mesh[i] += rhs.mesh[i];
  return *this;
}

RectMesh RectMesh::operator-() const
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  RectMesh temp(rows_y+2*nghost, cols_x+2*nghost);
  for (int i=0; i<len; i++)
    temp.mesh[i] = -1.0*mesh[i];
  return temp;
}

RectMesh RectMesh::operator-(const RectMesh& rhs) const
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);  
  RectMesh temp(rows_y +2*nghost, cols_x+2*nghost);
  for (int i=0; i<len; i++)
    temp.mesh[i] = mesh[i] - rhs.mesh[i];  
  return temp;
}

RectMesh& RectMesh::operator-=(const RectMesh& rhs) 
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  for (int i=0; i<len; i++)
    mesh[i] -= rhs.mesh[i];
  return *this;
}

RectMesh operator*(double scalar, const RectMesh& rhs)
{
  int rows = rhs.getrows() + 2*rhs.getnghost();
  int cols = rhs.getcols() + 2*rhs.getnghost();
  RectMesh temp(rows,cols);
  for (int i=0; i<rows*cols; i++)
    temp.mesh[i] = scalar*rhs.mesh[i];
  return temp;
}

RectMesh RectMesh::operator*(const RectMesh& rhs) const
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  RectMesh temp(rows_y+2*nghost, cols_x+2*nghost);
  for (int i=0; i<len; i++)
    temp.mesh[i] = mesh[i]*rhs.mesh[i];
  return temp;
}

RectMesh RectMesh::operator/(double scalar) const
{
  double val = 1./scalar;
  RectMesh temp(cols_x+2*nghost, rows_y+2*nghost);
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  for(int i=0; i<len; i++)
    temp.mesh[i] = mesh[i]*val;
  return temp;
}

/*----------------------- Other methods ------------------------*/

void RectMesh::fill(double value) const
{
  int rows = rows_y + 2*nghost;
  int cols = cols_x + 2*nghost;
  std::fill(mesh, mesh+rows*cols, value);
}

void RectMesh::print() const
{
  int row_end = rows_y + nghost;
  int col_end = cols_x + nghost;
  int jump = cols_x + 2*nghost;
  std::cout <<std::endl;
  for (int j = -nghost; j < row_end; j++)
    {
      for (int i = -nghost; i < col_end; i++)
	std::cout << std::setprecision(9)
		  << mesh [(i+nghost) + (j+nghost)*jump] << "\t";
      std::cout << std::endl;
    }
  std::cout <<std::endl;
}

void RectMesh::ln()
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  for(int i=0; i<len; i++)
    mesh[i] = log(mesh[i]); 
}

double RectMesh::sum() const
{
  int len = (cols_x+2*nghost)*(rows_y+2*nghost);
  double sum = 0.0;
  sum = std::accumulate(mesh, mesh+len, sum);
  return sum;
}

void RectMesh::writeH5(H5std_string filename) const
{
  hsize_t dims[2];
  const int RANK = 2;
  dims[0] = rows_y + 2*nghost; //slowest changing dimension
  dims[1] = cols_x + 2*nghost; //fastest changing dimension
  H5File file(filename, H5F_ACC_TRUNC);
  DataSpace dataspace(RANK,dims);
  const H5std_string dataset_name("RectMesh");
  DataSet dataset = file.createDataSet(dataset_name,
				       PredType::NATIVE_DOUBLE,
				       dataspace);
  dataset.write(mesh, PredType::NATIVE_DOUBLE);
}

void RectMesh::readH5(const char filename[])
{
  hsize_t DIM0 = cols_x + 2*nghost;
  hsize_t DIM1 = rows_y + 2*nghost;
  
  /* Open file and dataset using the default properties. */
  
  hid_t file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dset = H5Dopen (file, "RectMesh", H5P_DEFAULT);

  /* Check if the class instance has same dimensions with the dataset */

  hid_t dspace = H5Dget_space(dset);
  const int ndims = H5Sget_simple_extent_ndims(dspace);
  hsize_t dset_dims[ndims];
  H5Sget_simple_extent_dims(dspace, dset_dims, NULL);

  if (dset_dims[1] != DIM0 or dset_dims[0] != DIM1)
    {
      std::cout << "Dimensions do not agree. Cannot read from H5 file. Exiting."
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  /* If dimensions agree */
  
  herr_t status;
  H5D_layout_t layout;
  hsize_t dims[2] = {DIM0, DIM1};
  
  /* Retrieve the dataset creation property list, and print the 
     storage layout. */
  
  hid_t dcpl = H5Dget_create_plist (dset);
  layout = H5Pget_layout (dcpl);
  printf ("Storage layout for %s is: ", "RectMesh");
  switch (layout) {
  case H5D_COMPACT:
    printf ("H5D_COMPACT\n");
    break;
  case H5D_CONTIGUOUS:
    printf ("H5D_CONTIGUOUS\n");
    break;
  case H5D_CHUNKED:
    printf ("H5D_CHUNKED\n");
    break;
  case H5D_VIRTUAL:
    printf ("H5D_VIRTUAL\n");
    break;
  case H5D_LAYOUT_ERROR:
  case H5D_NLAYOUTS:
    printf ("H5D_LAYOUT_ERROR\n");
  }
  
  /* Read the data using the default properties. */
  
  status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    mesh);
  
  /*  Close and release resources. */
  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Fclose (file);
}
