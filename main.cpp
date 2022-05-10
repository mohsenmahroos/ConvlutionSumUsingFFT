// Accepted solution of FOUR ARRAYS problem 
// May 2022 Long-time Competitive Programming Contest 1
// at CodeChef.com

#include "FFT.h"
#include <iostream>
#include <numeric>

int main() {
    using namespace std;
    using i_limits = numeric_limits<size_t>;
    using ivector = vector<size_t>;
    using ull = unsigned long long;
    const size_t M = 4;
    ivector array_size(M);
    cin.tie(nullptr)->sync_with_stdio(false);
    for (size_t i = 0; i < M; ++i)
        cin >> array_size[i];
    ull K; cin >> K;
    const FFT_engine<18,ivector> FFT;
    const size_t N = FFT.N, u_min = i_limits::min(), u_max = i_limits::max();
    ivector f_min(M,u_max), f_max(M,u_min);
    vector<ivector> f(M,ivector(N));
    for (size_t i = 0; i < M; ++i)
        for (size_t elem; array_size[i] > 0; --array_size[i], ++f[i][elem])
            cin >> elem, 
            f_min[i] = min(f_min[i],elem),
            f_max[i] = max(f_max[i],elem);
    ivector g(N), h(N);
    FFT.convolution_sum(f[0],f[1],g);
    FFT.convolution_sum(f[2],f[3],h);
    for (size_t y = 1; y < N; ++y)
        h[y] += h[y-1];
    const auto y_index = [&](ull m, ull x) {
        const ull y = N-1;
        return x == 0 ? y: min(y,m/x); };
    const auto p_count = [&](ull m) {
        ull n = 0;
        for (size_t x = 0; x < N and n < K; ++x)
            if (g[x] > 0)
                n += 1ll*g[x]*h[y_index(m,x)];
        return n; };  
    const auto binary_search = [&](ull l, ull r) {
        ull ans = r;
        while (l <= r) {
            const ull m = (l+r)/2;
            if (p_count(m) < K)
                l = m+1;
            else
                r = m-1, ans = min(ans,m); }
        return ans; };
    const auto product_of_sums = [&](const ivector &x) {
        return 1ull*(x[0]+x[1])*(x[2]+x[3]); };
    cout << binary_search(product_of_sums(f_min),product_of_sums(f_max)); }
