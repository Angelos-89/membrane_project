#include "Site.hpp"

bool operator==(const Site& lhs, const Site& rhs)
{
  return ( (lhs.getx() == rhs.getx()) and (lhs.gety() == rhs.gety()) );
}

