#include"polynomial.h"
using namespace std;

template<typename Z> pair<polynomial, polynomial> polynomial::division(polynomial const &b) const{
   assert(m == b.m);
   int d1 = deg(), d2 = b.deg();
   assert(d2 >= 0);
   if(d1 < d2){
      auto r = *this;
      return {polynomial(m), r.trunc()};
   }
   if(d2 <= 32){
      return osoi_warizan(b);
   }
   auto q = *this/b, r = *this-b*q;
   return {q, r};
}
