/* This class represents a 3D mathematical vector. 
   Its coefficients are of type double. */

#ifndef VEC3D_HPP
#define VEC3D_HPP

#include <iostream>

class Vec3d
{
private:
  
  double x,y,z;

public:

  /* Constructors */
  Vec3d(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {};
  
  /* Algebra */
  Vec3d operator+(const Vec3d& aVector) const
  {
    return Vec3d{x + aVector.x, y + aVector.y, z + aVector.z};    
  }
  Vec3d operator-(const Vec3d& aVector) const
  {
    return Vec3d{x - aVector.x, y - aVector.y, z - aVector.z};
  }
  Vec3d operator*(double scalar) const
  {
    return Vec3d{x*scalar, y*scalar, z*scalar};
  }
  friend Vec3d operator*(double scalar, Vec3d& aVector) 
  {
    return Vec3d{scalar*aVector.x, scalar*aVector.y, scalar*aVector.z};
  }
  Vec3d operator/(double scalar) const
  {
    return Vec3d{x/scalar, y/scalar, z/scalar};
  }

  /* Othe methods */
  double X() const {return x;}
  double Y() const {return y;}
  double Z() const {return z;}
  void print() const
  {
    std::cout << "(" << x << "," << y << "," << z << ")" <<std::endl;
  }
};
#endif 
