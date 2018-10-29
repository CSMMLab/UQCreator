#ifndef VECTOR_H
#define VECTOR_H

#include <sstream>
#include <cmath>

namespace VectorSpace {

template<class T>
class Vector
{
private:
   unsigned _N;
   T* _data;

public:
   typedef T* iterator;
   typedef const T* const_iterator;

   Vector();
   Vector(unsigned n, bool skipZeroInit = false);
   Vector(unsigned n, T init);
   Vector(const Vector& other);
   ~Vector();

   T& operator[](unsigned i);
   const T& operator[](unsigned i) const;
   Vector operator+(const Vector& other) const;
   Vector operator-(const Vector& other) const;
   Vector operator*(const Vector& other) const;
   Vector operator*(const double& scalar) const;
   Vector operator/(const Vector& other) const;
   void operator=(const Vector& other) const;

   unsigned size() const;
   void reset();
   void resize(unsigned newSize);
   iterator begin();
   iterator end();
};

}

#endif // VECTOR_H

