/* This header file contains all the functions that 
   have objects of the class Vec3d as argument. */

#ifndef VEC3DLIB_HPP
#define VEC3DLIB_HPP

#include "Vec3d.hpp"

double Dot(const Vec3d& v, const Vec3d& w);          
Vec3d Cross(const Vec3d& v, const Vec3d& w);      
double Square(const Vec3d& v);                          
double Norm(const Vec3d& v);                            
Vec3d Normalize(const Vec3d& v, double magnitude);   
double Distance(const Vec3d& v, const Vec3d& w);

#endif 
