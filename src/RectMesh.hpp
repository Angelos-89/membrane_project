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

  void CheckConstr(int& cols, int& rows, int& ghost);
  RectMesh(int cols=2, int rows=2, int ghost=0);

  /* Rule of 5 */
  ~RectMesh();
  RectMesh(const RectMesh& obj);
  RectMesh(RectMesh&& tmp) noexcept;
  RectMesh& operator=(const RectMesh& rhs);
  RectMesh& operator=(RectMesh&& tmp);

  /* Operator () overloading */
  void CheckIndex(int& i, int& j) const;
  double& operator()(int i, int j) const;
  double& operator()(int i, int j);

  /* Other methods */
  void print() const;
  int getcols() const {return cols_x;}
  int getrows() const {return rows_y;}
  int getnghost() const {return nghost;}
  double* getmemory() const {return mesh;}
  void ln();
  double sum() const;
  void write(H5std_string filename);
  
};
#endif
