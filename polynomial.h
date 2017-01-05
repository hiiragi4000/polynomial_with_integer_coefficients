#include<cmath>
#include<cassert>
#include<complex>
#include<vector>
#include<type_traits>

using namespace std;

#define FFT_LEN 131072
#ifndef M_PI
const double M_PI = acos(-1);
#endif

vector<complex<double>> fft(const vector<complex<double>> &x, const int &inv){
    static bool fft_ready = false;
    static complex<double> loli[FFT_LEN];
    int n = x.size();
    assert((n&-n)==n && n<=FFT_LEN && abs(inv)==1);
    if(!fft_ready){
        for(int k=0; k<FFT_LEN; k++){
            loli[k] = exp(complex<double>(0, 2*M_PI*k/FFT_LEN));
        }
        fft_ready = true;
    }
    vector<complex<double>> X = x;
    for(int i=1, j=0; i<n; i++){
        for(int k=n>>1; !((j^=k)&k); k>>=1);
        if(i < j){
            swap(X[i], X[j]);
        }
    }
    for(int i=2; i<=n; i*=2){
        int d = (inv==1)? FFT_LEN-(FFT_LEN/i): FFT_LEN/i;
        for(int j=0; j<n; j+=i){
            for(int k=0, a=0; k<i/2; k++, a=(a+d)%FFT_LEN){
                complex<double> s = X[j+k], t = loli[a] * X[j+k+i/2];
                X[j+k] = s + t;
                X[j+k+i/2] = s - t;
            }
        }
    }
    if(inv == -1){
        for(int i=0; i<(int)X.size(); i++){
            X[i] /= n;
        }
    }
    return X;
}

template<class R> class polynomial{
private:
    vector<R> a;
    polynomial<R> karatsuba(const polynomial<R> &b) const{
        if(!size() || !b.size()){
            return polynomial<R>();
        }
        polynomial<R> result(size()+b.size()-1);
        if((int)size()<=64 || (int)b.size()<=64){
            for(int i=0; i<(int)size(); i++) for(int j=0; j<(int)b.size(); j++){
                result[i+j] += a[i] * b[j];
            }
        }else{
            int n = max(size(), b.size());
            polynomial<R> alpha(max((int)size()-n/2, 0)), beta(min((int)size(), n/2)), gamma(max((int)b.size()-n/2, 0)), delta(min((int)b.size(), n/2));
            for(int i=0; i<(int)alpha.size(); i++){
                alpha[i] = a[n/2+i];
            }
            for(int i=0; i<(int)beta.size(); i++){
                beta[i] = a[i];
            }
            for(int i=0; i<(int)gamma.size(); i++){
                gamma[i] = b[n/2+i];
            }
            for(int i=0; i<(int)delta.size(); i++){
                delta[i] = b[i];
            }
            polynomial<R> p1 = alpha*gamma, p2 = beta*delta, p3 = (alpha+beta)*(gamma+delta) - p1 - p2;
            for(int i=0; i<(int)p1.size(); i++){
                result[(n-n%2)+i] += p1[i];
            }
            for(int i=0; i<(int)p2.size(); i++){
                result[i] += p2[i];
            }
            for(int i=0; i<(int)p3.size(); i++){
                result[n/2+i] += p3[i];
            }
        }
        return result;
    }
    template<class S> class enable_if<is_same<R, S>::value && is_same<S, complex<double>>::value, polynomial<S>>::type kakeru(const polynomial<S> &b) const{
        int n = size()+b.size()-1;
        if((int)size()<=16 || (int)b.size()<=16){
            return karatsuba(b);
        }
        for(; (n&-n)!=n; n+=n&-n);
        vector<complex<double>> x(n), y(n);
        for(int i=0; i<(int)size(); i++){
            x[i] = a[i];
        }
        for(int i=0; i<(int)b.size(); i++){
            y[i] = b[i];
        }
        x = fft(x, 1);
        y = fft(y, 1);
        for(int i=0; i<n; i++){
            x[i] *= y[i];
        }
        polynomial<complex<double>> result(fft(x, -1));
        result.resize(size()+b.size()-1);
        return result;
    }
    template<class S> class enable_if<is_same<R, S>::value && is_arithmetic<S>::value, polynomial<S>>::type kakeru(const polynomial<S> &b) const{
        polynomial<complex<double>> f(size()), g(b.size());
        for(int i=0; i<(int)size(); i++){
            f[i] = a[i];
        }
        for(int i=0; i<(int)b.size(); i++){
            g[i] = b[i];
        }
        polynomial<complex<double>> h = f * g;
        polynomial<S> result(h.size());
        if(is_floating_point<S>::value){
            for(int i=0; i<(int)h.size(); i++){
                result[i] = h[i].real();
            }
        }else{
            for(int i=0; i<(int)h.size(); i++){
                result[i] = (S)floor(h[i].real()+0.5);
            }
        }
        return result;
    }
    template<class S> class enable_if<is_same<R, S>::value && !is_same<S, complex<double>>::value && !is_arithmetic<S>::value>::type kakeru(const polynomial<S> &b) const{
        return karatsuba(b);
    }
public:
    polynomial(const size_t &n = 0){
        a = vector<R>(n);
    }
    polynomial(const vector<R> &coef){
        a = coef;
    }
    size_t size() const{
        return a.size();
    }
    void resize(const size_t &n){
        a.resize(n);
    }
    R& operator [](const int &i){
        assert(0<=i && i<(int)a.size());
        return a[i];
    }
    const R& operator [](const int &i) const{
        assert(0<=i && i<(int)a.size());
        return a[i];
    }
    polynomial<R> operator +(const polynomial<R> &b) const{
        polynomial<R> result = *this;
        if(size() < b.size()){
            result.resize(b.size());
        }
        for(int i=0; i<(int)b.size(); i++){
            result[i] += b[i];
        }
        return result;
    }
    polynomial<R> operator +=(const polynomial<R> &b){
        return (*this) = (*this) + b;
    }
    polynomial<R> operator -() const{
        polynomial<R> result = *this;
        for(int i=0; i<(int)size(); i++){
            result[i] = -result[i];
        }
        return result;
    }
    polynomial<R> operator -(const polynomial<R> &b) const{
        return (*this) + (-b);
    }
    polynomial<R> operator -=(const polynomial<R> &b){
        return (*this) = (*this) - b;
    }
    polynomial<R> operator *(const R &c) const{
        polynomial<R> result = (*this);
        for(int i=0; i<(int)result.size(); i++){
            result[i] *= c;
        }
        return result;
    }
    polynomial<R> operator *=(const R &c){
        return (*this) = (*this) * c;
    }
    polynomial<R> operator *(const polynomial<R> &b) const{
        return kakeru(b);
    }
    polynomial<R> operator *=(const polynomial<R> &b){
        return (*this) = (*this) * b;
    }
    polynomial<R> operator /(const R &c) const{
        polynomial<R> result = (*this);
        for(int i=0; i<(int)result.size(); i++){
            result[i] /= c;
        }
        return result;
    }
    polynomial<R> operator /=(const R &c){
        return (*this) = (*this) / c;
    }
    polynomial<R> operator /(const polynomial<R> &b) const{
        return divide(b).first;
    }
    polynomial<R> operator /=(const polynomial<R> &b){
        return (*this) = (*this) / b;
    }
    polynomial<R> operator %(const polynomial<R> &b) const{
        return divide(b).second;
    }
    polynomial<R> operator %=(const polynomial<R> &b){
        return (*this) = (*this) % b;
    }
    pair<polynomial<R>, polynomial<R>> divide(const polynomial<R> &b) const{
        assert(b.size());
        if(b.size() > size()){
            return make_pair(polynomial<R>(), *this);
        }
        polynomial<R> q(size()-b.size()+1), r = (*this);
        for(int i=q.size()-1; i>=0; i--){
            q[i] = r[b.size()+i-1] / b[b.size()-1];
            for(int j=0; j<(int)b.size(); j++){
                r[i+j] -= q[i] * b[j];
            }
        }
        r.resize(b.size()-1);
        return make_pair(q, r);
    }
};
template<class R> polynomial<R> operator *(const R &c, polynomial<R> f){
    return f*c;
}
