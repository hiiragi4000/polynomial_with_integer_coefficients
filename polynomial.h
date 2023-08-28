#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include<algorithm>
#include<complex>
#include<type_traits>
#include<utility>
#include<vector>
#include<cassert>
#include<cmath>

// This function is NOT atomic. Do NOT use it in any concurrent situation.
inline std::vector<std::complex<double>> fft(std::vector<std::complex<double>> const &x, int inv){
   constexpr double PI = 3.14159265358979323846;
   static int fft_len = 0;
   static std::vector<std::complex<double>> table;
   size_t n = x.size();
   assert((n&~n+1)==n && std::abs(inv)==1);
   if(fft_len < n){
      table.resize(fft_len = (int)n);
      for(size_t k=0; k<=n-1; ++k){
         table[k] = std::exp(std::complex<double>(0, 2*k*PI/n));
      }
   }
   auto X = x;
   for(size_t i=1, j=0; i<=n-1; ++i){
      for(size_t k=n>>1; !((j^=k)&k); k>>=1);
      if(i < j){
         std::swap(X[i], X[j]);
      }
   }
   for(size_t i=2; i<=n; i*=2){
      int d = (inv==1)? fft_len-(fft_len/(int)i): fft_len/(int)i;
      for(size_t j=0; j<n; j+=i){
         for(size_t k=0, a=0; k<i/2; ++k, a=(a+d)%fft_len){
            auto s = X[j+k], t = table[a]*X[j+k+i/2];
            X[j+k] = s+t;
            X[j+k+i/2] = s-t;
         }
      }
   }
   if(inv == -1){
      for(auto &Xi: X){
         Xi /= (int)n;
      }
   }
   return X;
}

#define I64 long long

// This function is NOT atomic. Do NOT use it in any concurrent situation.
#define NTT_LEN (1<<27)
constexpr I64 nttp = (15<<27)+1;
inline std::vector<I64> ntt(std::vector<I64> const &x, int inv){
   constexpr I64 pr = 31;
   static int ntt_len = 0;
   static std::vector<I64> table;
   size_t n = x.size();
   assert((n&~n+1)==n && n<=NTT_LEN && std::abs(inv)==1);
   for(size_t i=0; i<=n-1; ++i){
      assert(0<=x[i] && x[i]<=nttp-1);
   }
   if(ntt_len < n){
      table.resize(ntt_len = (int)n);
      table[0] = 1;
      if(n >= 2){
         table[1] = 1;
         for(int T=15; T-->0; ){
            table[1] = pr*table[1]%nttp;
         }
         for(int i=NTT_LEN/(int)n; i>1; i>>=1){
            table[1] = table[1]*table[1]%nttp;
         }
         for(size_t i=2; i<=n-1; ++i){
            table[i] = table[1]*table[i-1]%nttp;
         }
      }
   }
   auto X = x;
   for(size_t i=1, j=0; i<=n-1; ++i){
      for(size_t k=n>>1; !((j^=k)&k); k>>=1);
      if(i < j){
         std::swap(X[i], X[j]);
      }
   }
   for(size_t i=2; i<=n; i*=2){
      int d = (inv==1)? ntt_len-(ntt_len/(int)i): ntt_len/(int)i;
      for(size_t j=0; j<n; j+=i){
         for(size_t k=0, a=0; k<i/2; ++k, a=(a+d)%ntt_len){
            auto s = X[j+k], t = table[a]*X[j+k+i/2]%nttp;
            X[j+k] = (s+t)%nttp;
            X[j+k+i/2] = (s+nttp-t)%nttp;
         }
      }
   }
   if(inv == -1){
      for(auto &Xi: X){
         Xi = Xi*(nttp-(nttp-1)/n)%nttp;
      }
   }
   return X;
}
#undef NTT_LEN

template<typename Z> struct polynomial{
   static_assert(std::is_integral<Z>::value, "Usage: polynomial<Z> for an integer type Z.\n");
   explicit polynomial(Z const &m=0) noexcept: m(m){
      assert(m>=0 && m!=1);
   }
   polynomial(std::vector<Z> coef, Z const &m=0) noexcept: a(std::move(coef)), m(m){
      assert(m>=0 && m!=1);
   }
   void resize(size_t n){
      a.resize(n);
   }
   size_t size() const noexcept{
      return a.size();
   }
   Z &operator[](size_t i){
      if(i >= a.size()) a.resize(i+1);
      return a[i];
   }
   Z operator[](size_t i) const noexcept{
      return i<a.size()? a[i]: 0;
   }
   Z mod() const noexcept{
      return m;
   }
   int deg() const noexcept{
      int d = (int)size()-1;
      if(m){
         for(; d>=0&&a[d]%m==0; --d);
      }else{
         for(; d>=0&&!a[d]; --d);
      }
      return d;
   }
   polynomial &trunc(){
      if(m){
         for(auto &ai: a){
            ai = (ai%m+m)%m;
         }
      }
      while(!a.empty() && !a.back()){
         a.pop_back();
      }
      return *this;
   }
   polynomial operator+() const noexcept{
      return *this;
   }
   polynomial operator-() const noexcept{
      auto res = *this;
      for(auto &ri: res.a){
         ri = -ri;
      }
      return res.trunc();
   }
   polynomial &operator+=(polynomial const &b){
      assert(m == b.m);
      if(size() < b.size()){
         resize(b.size());
      }
      for(size_t i=0; i<size(); ++i){
         a[i] += b[i];
      }
      return trunc();
   }
   polynomial &operator-=(polynomial const &b){
      return *this += -b;
   }
   polynomial &operator*=(Z const &c){
      for(auto &ai: a){
         ai *= c;
      }
      return trunc();
   }
   polynomial &operator*=(polynomial const &b){
      assert(m == b.m);
      auto d1 = deg(), d2 = b.deg();
      if(d1<=16 || d2<=16){
         return *this = osoi_kakezan(b);
      }
      size_t n = d1+d2+1;
      resize(n);
      for(; (n&~n+1)!=n; n+=n&~n+1);
      std::vector<std::complex<double>> x(n), y(n);
      for(int i=0; i<=d1; ++i){
         x[i] = a[i];
      }
      for(int i=0; i<=d2; ++i){
         y[i] = b[i];
      }
      x = fft(x, 1); y = fft(y, 1);
      for(size_t i=0; i<x.size(); ++i){
         x[i] *= y[i];
      }
      x = fft(x, -1);
      std::vector<I64> xi(n), eta(n);
      for(int i=0; i<=d1; ++i){
         xi[i] = (a[i]%nttp+nttp)%nttp;
      }
      for(int i=0; i<=d2; ++i){
         eta[i] = (b[i]%nttp+nttp)%nttp;
      }
      xi = ntt(xi, 1); eta = ntt(eta, 1);
      for(size_t i=0; i<xi.size(); ++i){
         xi[i] = xi[i]*eta[i]%nttp;
      }
      xi = ntt(xi, -1);
      for(int i=0; i<=d1+d2; ++i){
         I64 q = (I64)floor((x[i].real()-xi[i])/nttp+.5);
         if(!m){
            a[i] = (Z)(q*nttp+xi[i]);
         }else{
            a[i] = (Z)((q%m)*(nttp%m)%m+xi[i])%m;
         }
      }
      return trunc();
   }
   polynomial &operator/=(Z const &c){
      for(size_t i=0; i<size(); ++i){
         a[i] /= c;
      }
      return *this;
   }
   polynomial &operator%=(Z const &c){
      for(size_t i=0; i<size(); ++i){
         a[i] %= c;
      }
      return *this;
   }
   polynomial rev() const{
      auto res = *this;
      res.trunc();
      reverse(res.a.begin(), res.a.end());
      return res.trunc();
   }
   polynomial invmodxn(int n) const{
      assert(n>0 && deg()>=0 && abs(a[0])==1);
      if(n == 1){
         polynomial b(m);
         b[0] = a[0]; return b;
      }
      int l = (n+1)/2;
      polynomial p = invmodxn(l), a0(m), a1(m);
      a0.resize(l); a1.resize(l);
      for(int i=0; i<=l-1 && i<(int)size(); ++i){
         a0[i] = a[i];
      }
      a0.trunc();
      for(int i=0; i<=l-1 && l+i<(int)size(); ++i){
         a1[i] = a[l+i];
      }
      a1.trunc();
      polynomial a0p = a0*p, b(m);
      b.resize(l);
      for(int i=0; i<=l-1 && l+i<(int)a0p.size(); ++i){
         b[i] = a0p[l+i];
      }
      b.trunc();
      auto q = a1*p+b;
      q.resize(l); q *= -p;
      polynomial res(m); res.resize(n);
      for(int i=0; i<=l-1 && i<(int)p.size(); ++i){
         res[i] = p[i];
      }
      for(int i=l; i<=n-1 && i-l<(int)q.size(); ++i){
         res[i] = q[i-l];
      }
      return res.trunc();
   }
   polynomial &operator/=(polynomial const &b){
      assert(m == b.m);
      int d1 = deg(), d2 = b.deg();
      assert(d2 >= 0);
      if(d1 < d2){
         return *this = polynomial{m};
      }
      if(d2 <= 32){
         return *this = osoi_warizan(b).first;
      }
      auto rf = rev(), rg = b.rev();
      *this = rf*rg.invmodxn(d1-d2+1);
      resize(d1-d2+1);
      std::reverse(a.begin(), a.end());
      return *this;
   }
   polynomial &operator%=(polynomial const &b){
      assert(m == b.m);
      int d1 = deg(), d2 = b.deg();
      assert(d2 >= 0);
      if(d1 < d2){
         return *this;
      }
      if(d2 <= 32){
         return *this = osoi_warizan(b).second;
      }
      return *this = *this - ((*this)/b)*b;
   }
   std::pair<polynomial, polynomial> division(polynomial const &b) const;
private:
   std::vector<Z> a;
   Z m;
   polynomial osoi_kakezan(polynomial const &b) const{
      assert(m == b.m);
      int d1 = deg(), d2 = b.deg();
      if(d1<0 || d2<0){
         return polynomial(m);
      }
      polynomial res(m); res.resize(d1+d2+1);
      for(int i=0; i<=d1; ++i) for(int j=0; j<=d2; ++j){
         if(m){
            res[i+j] = (res[i+j]+a[i]*b[j]%m)%m;
         }else{
            res[i+j] += a[i]*b[j];
         }
      }
      return res.trunc();
   }
   std::pair<polynomial, polynomial> osoi_warizan(polynomial const &b) const{
      assert(m == b.m);
      auto r = *this, g = b;
      r.trunc(); g.trunc();
      int d1 = r.deg(), d2 = g.deg();
      assert(d2>=0 && std::abs(g[d2])==1);
      if(d1 < d2){
         return {polynomial(m), r};
      }
      polynomial q(m); q.resize(d1-d2+1);
      for(int i=d1-d2; i>=0; --i){
         q[i] = g[d2]*r[d2+i];
         for(int j=0; j<=d2; ++j){
            if(m){
               r[i+j] = (r[i+j]+(m-q[i])*g[j]%m)%m;
            }else{
               r[i+j] -= q[i]*g[j];
            }
         }
      }
      return {q, r.trunc()};
   }
};
template<typename Z> std::pair<polynomial<Z>, polynomial<Z>> polynomial<Z>::division(polynomial<Z> const &b) const{
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
template<typename Z> polynomial<Z> operator+(polynomial<Z> f, polynomial<Z> const &g){
   return f += g;
}
template<typename Z> polynomial<Z> operator-(polynomial<Z> f, polynomial<Z> const &g){
   return f -= g;
}
template<typename Z> polynomial<Z> operator*(Z const &c, polynomial<Z> f){
   return f *= c;
}
template<typename Z> polynomial<Z> operator*(polynomial<Z> f, Z const &c){
   return f *= c;
}
template<typename Z> polynomial<Z> operator*(polynomial<Z> f, polynomial<Z> const &g){
   return f *= g;
}
template<typename Z> polynomial<Z> operator/(polynomial<Z> f, Z const &c){
   return f /= c;
}
template<typename Z> polynomial<Z> operator/(polynomial<Z> f, polynomial<Z> const &g){
   return f /= g;
}
template<typename Z> polynomial<Z> operator%(polynomial<Z> f, Z const &c){
   return f %= c;
}
template<typename Z> polynomial<Z> operator%(polynomial<Z> f, polynomial<Z> const &g){
   return f %= g;
}
template<typename Z> bool operator==(polynomial<Z> const &f, polynomial<Z> const &g){
   if(f.mod()!=g.mod() || f.deg()!=g.deg()){
      return false;
   }
   auto m = f.mod();
   int n = f.deg();
   if(m == 0){
      for(int i=0; i<=n; ++i){
         if(f[i] != g[i]){
            return false;
         }
      }
   }else{
      for(int i=0; i<=n; ++i){
         if(f[i]%m != g[i]%m){
            return false;
         }
      }
   }
   return true;
}
template<typename Z> bool operator!=(polynomial<Z> const &f, polynomial<Z> const &g){
   return !(f == g);
}

#undef I64

#endif
