#include"polynomial.h"
#include<cstdio>
using namespace std;

bool test_small_01(){
   // f(x) = x^2+2x+2, g(x) = -x^2+2x-2
   polynomial<int> f({2, 2, 1}), g({-2, 2, -1});
   bool res = true;
   // (f+g)(x) = 4x
   res &= f+g == polynomial<int>({0, 4, 0, 0, 0});
   // (f-g)(x) = 2x^2+4
   res &= f-g == polynomial<int>({4, 0, 2, 0, 0});
   // (fg)(x) = -x^4-4
   res &= f*g == polynomial<int>({-4, 0, 0, 0, -1});
   // f(x) = (-1)g(x) + (4x)
   res &= f/g == polynomial<int>({-1, 0, 0, 0, 0});
   res &= f%g == polynomial<int>({0, 4, 0, 0, 0});
   return res;
}

bool test_large_01(){
   // f(x) = x^20000+x^10000+1, g(x) = x^25000-x^15000+x^5000
   polynomial<int> f, g;
   f[20000] = f[10000] = f[0] = 1;
   g[25000] = g[5000] = 1, g[15000] = -1;
   bool res = true;
   // (f+g)(x) = x^25k + x^20k - x^15k + x^10k + x^5k + 1
   polynomial<int> h;
   h[25000] = h[20000] = h[10000] = h[5000] = h[0] = 1;
   h[15000] = -1;
   res &= f+g == h;
   // (f-g)(x) = -x^25k + x^20k + x^15k + x^10k - x^5k + 1
   h[25000] = h[5000] = -1, h[15000] = 1;
   res &= f-g == h;
   // (fg)(x) = x^45k + x^25k + x^5k
   polynomial<int> fg;
   fg[45000] = fg[25000] = fg[5000] = 1;
   res &= f*g == fg;
   // f(x) = (0)g(x) + f(x)
   polynomial<int> zero;
   res &= f/g == zero;
   res &= f%g == f;
   // g(x) = (x^5k)f(x) + (-2x^15k)
   polynomial<int> x5k, x15k;
   x5k[5000] = 1, x15k[15000] = 1;
   res &= g/f == x5k;
   res &= g%f == -2*x15k;
   return res;
}

int main(){
   bool res = true;
   res &= test_small_01();
   res &= test_large_01();
   if(res){
      puts("All tests passed.");
   }else{
      puts("Some tests failed.");
      return 1;
   }
}
