/* A class representing a 2D scalar field defined on a 
   rectangular uniform grid (mesh or lattice). */

#ifndef RECTMESH_HPP
#define RECTMESH_HPP

#include <fftw3.h>
#include "H5Cpp.h"
using namespace H5;

class RectMesh
{
private:
  
  int cols_x;
  int rows_y;
  int nghost;
  double* mesh;

public:

  /* Constructor */
  void CheckConstr(int& cols, int& rows, int& ghost);
  explicit RectMesh(int cols=2, int rows=2, int ghost=0);

  /* Rule of 5 */
  ~RectMesh();
  RectMesh(const RectMesh& obj);
  RectMesh(RectMesh&& tmp) noexcept;
  RectMesh& operator=(const RectMesh& rhs);
  RectMesh& operator=(RectMesh&& tmp);

  /* Operator () overloading */
  double& operator()(int i, int j) const;
  double& operator()(int i, int j);

  /* Algebra */
  RectMesh operator+(const RectMesh& rhs) const;
  RectMesh& operator+=(const RectMesh& rhs);
  RectMesh operator-() const;
  RectMesh operator-(const RectMesh& rhs) const;
  RectMesh& operator-=(const RectMesh& rhs);
  friend RectMesh operator*(double scalar, const RectMesh& rhs);
  RectMesh operator*(const RectMesh& rhs) const;
  RectMesh operator/(double scalar) const;

  /* Other methods */
  void fill(double value) const;
  void print() const;
  int getcols() const {return cols_x;}
  int getrows() const {return rows_y;}
  int getnghost() const {return nghost;}
  double* getmemory() const {return mesh;}
  void ln();
  double sum() const;
  void writeH5(H5std_string filename) const;
  void readH5(const char filename[]);  
};
#endif
