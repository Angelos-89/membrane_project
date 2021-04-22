/* Class representing a lattice site (i,j). Here i is an int that indicates
   the x coordinate and j indicates the y coordinate. */

#ifndef SITE_HPP
#define SITE_HPP

#include <iostream>

class Site
{
private:

  int x;
  int y;

public:

  Site(int xx=0, int yy=0) : x(xx), y(yy) {};
  void set(int xx, int yy){ x = xx; y=yy;}
  friend bool operator==(const Site& lhs, const Site& rhs);
  int getx() const {return x;}
  int gety() const {return y;}
  void print() const
  {
    std::cout << std::endl;
    std::cout << "(" << x << "," << y << ")" << std::endl;
    std::cout << std::endl;
  }

  double dist(const Site& a_site);
  
};
#endif
