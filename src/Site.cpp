#include "Site.hpp"
#include <math.h>
#include <iostream>

bool operator==(const Site& lhs, const Site& rhs)
{
  return ( (lhs.getx() == rhs.getx()) and (lhs.gety() == rhs.gety()) );
}

double Site::dist(const Site& a_site)
{
  int x1 = (*this).getx();
  int y1 = (*this).gety();
  int x2 = a_site.getx();
  int y2 = a_site.gety();
  return sqrt( pow(x2-x1,2) + pow(y2-y1,2) );
}
