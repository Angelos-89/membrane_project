#include <iostream>
#include <cmath>
#include "Vec3d.hpp"
#include "Vec3dLib.hpp"

double Dot(const Vec3d& v, const Vec3d& w) 
{
  return v.X()*w.X() + v.Y()*w.Y() + v.Z()*w.Z();
}

Vec3d Cross(const Vec3d& v, const Vec3d& w) 
{
  return {v.Y()*w.Z()-v.Z()*w.Y(),
	  v.Z()*w.X()-v.X()*w.Z(),
	  v.X()*w.Y()-v.Y()*w.X()};
}

double Square(const Vec3d& vec1) 
{
  return Dot(vec1,vec1);
}

double Norm(const Vec3d& vec1)
{
  return sqrt(Dot(vec1,vec1));
}

Vec3d Normalize(const Vec3d& v, double magnitude)
{
  if (magnitude == 0)
    {
      std::cout << "Division with zero! Vector cannot be normalized."
		<< std::endl; 
      return v;
    }
  else
    {
      return {v.X()/magnitude, v.Y()/magnitude, v.Z()/magnitude};
    }
}

double Distance(Vec3d& v, Vec3d& w)
{
  Vec3d diff = v-w;
  return Norm(diff);
}
