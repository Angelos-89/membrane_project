#include "Site.hpp"

bool operator==(const Site& lhs, const Site& rhs)
{
  bool var;
  if (lhs.getx() == rhs.getx() && lhs.gety() == rhs.gety())
    return (bool) 1;
  else
    return (bool) 0;
}

