#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "RectMesh.hpp"

/*--------------------- Constructor -----------------------------*/
void RectMesh::CheckConstr(int& cols, int& rows, int& ghost)
{
  if (cols < 2 || rows < 2 || ghost < 0)
    {
      std::cout << "The smallest mesh must have dimensions (2x2)"
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

/*------- Copy assignment operator -------*/
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
void RectMesh::CheckIndex(int& i, int& j) const
{
  if ( i < -nghost || i >= cols_x+nghost ||
       j < -nghost || j >= rows_y+nghost )
    {
      std::cout << "Index out of bounds. Exiting."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
}

double& RectMesh::operator()(int i, int j) const
{
  CheckIndex(i,j);
  int jump = cols_x + 2*nghost;
  return mesh[ (i+nghost) + (j+nghost)*jump ];
}

double& RectMesh::operator()(int i, int j)
{
  CheckIndex(i,j);
  int jump = cols_x + 2*nghost;
  return mesh[ (i+nghost) + (j+nghost)*jump ];
}

/*------------------------- Other methods ---------------------*/
void RectMesh::print() const
{
  int row_end = rows_y + nghost;
  int col_end = cols_x + nghost;
  int jump = cols_x + 2*nghost;
  std::cout <<std::endl;
  for (int j = -nghost; j < row_end; j++)
    {
      for (int i = -nghost; i < col_end; i++)
	std::cout << mesh [(i+nghost) + (j+nghost)*jump] << "\t";
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

void RectMesh::write(H5std_string filename)
{
  hsize_t dims[2];
  const int RANK = 2;
  dims[0] = cols_x + 2*nghost;
  dims[1] = rows_y + 2*nghost;
  H5File file(filename, H5F_ACC_TRUNC);
  DataSpace dataspace(RANK,dims);
  const H5std_string dataset_name("RectMesh");
  DataSet dataset = file.createDataSet(dataset_name,
				       PredType::NATIVE_DOUBLE,
				       dataspace);
  dataset.write(mesh, PredType::NATIVE_DOUBLE);
}
