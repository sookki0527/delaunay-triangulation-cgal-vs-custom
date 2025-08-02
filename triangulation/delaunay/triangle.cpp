#include "triangle.h"

namespace krs{
    
  template<typename U>
  std::ostream&
      operator <<(std::ostream& str, const Triangle<U>& t)
  {
      return str << "Triangle:" << "\n\t" <<
          *t.a << "\n\t" <<
          *t.b << "\n\t" <<
          *t.c << '\n';
  }
  
}